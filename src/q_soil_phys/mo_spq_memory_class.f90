!> QUINCY soil-physics process memory
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of (memory) variables for the soil-physics-quincy process
!>
MODULE mo_spq_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE, LAKE_TYPE, BARE_TYPE, GLACIER_TYPE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_spq_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars      = 240


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for spq memory
  !!
  !! @par includes: \n
  !!    spq variables spq\%VAR
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_spq_memory


    TYPE(t_jsb_var_real2d)          :: &
                                    soil_depth      , & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
                                    root_depth      , & !< actual rooting depth [m]
                                    t_srf_new       , & !< current surface temperature [K]
                                    t_srf_old       , & !< previous surface temperature [K]
                                    temp_srf_eff_4  , & !< effective surface temperature ** 4.0 [K**4.0] (in jsbach: t_eff4)
                                    zril_old        , & !< previous timestep's Reynold's number
                                    elevation       , & !< currently 0.0 [m]
                                    qsat_star       , & !< saturation specific humidity at surface [kg kg-1 ?]
                                    s_star              !< surface dry static energy               [m2 s-2 ?]

    TYPE(t_jsb_var_real2d)          :: &
      & evapotranspiration                !< evapotranspiration [kg m-2 s-1]

    TYPE(t_jsb_var_real2d)          :: &
      &   num_sl_above_bedrock                          !< number of soil layers above bedrock, i.e., with a width larger eps8

    TYPE(t_jsb_var_real2d)          :: &
      & q_sat_srf                          !< surface saturated specific humidity  [?]

    ! JSBACH SOIL
    TYPE(t_jsb_var_real3d)          :: &
                                    soil_depth_sl, &            !< Width ! of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_width_sl, &        !< width of each soil layer that can be saturated with water, identical with soil_depth_sl (until bedrock) [m]
                                    soil_lay_depth_center_sl, & !< depth at the center of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_depth_ubound_sl, & !< depth at the upper bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]
                                    soil_lay_depth_lbound_sl    !< depth at the lower bounddary (ubound > lbound !) of each soil layer that can be saturated with water (until bedrock) [m]

    ! JSBACH HYDRO MEMORY (temporary here to get rid of vegsoil construct)
    TYPE(t_jsb_var_real2d)          :: &
                                    interception,     & !< surface interception evaporation                                         [kg m-2 s-1]
                                    evapopot,         & !< potential evaporation                                                    [kg m-2 s-1]
                                    evaporation,      & !< surface evaporation                                                      [kg m-2 s-1]
                                    evaporation_snow, & !< snow surface evaporation                                                 [kg m-2 s-1] unit correct ?
                                    srf_runoff,       & !< surface runoff                                                           [kg m-2 s-1]
                                    drainage,         & !< drainage                                                                 [kg m-2 s-1]
                                    gw_runoff,        & !< lateral runoff                                                           [kg m-2 s-1]
                                    drainage_fraction   !< water mass fraction lost from the soil column by drainage                [1 s-1]

    ! Additional variables for spq
    TYPE(t_jsb_var_real2d)          :: &
                                    w_skin,           & !< Water content in skin reservoir                                          [m]
                                    w_soil_root,      & !< Water content in root zone of the soil                                   [m]
                                    w_soil_root_fc,   & !< Water content at field capacity in root zone of the soil                 [m]
                                    w_soil_root_pwp     !< Water content at permanent wilting point in root zone of the soil        [m]

    TYPE(t_jsb_var_real3d)          :: &
      & frac_w_lat_loss_sl, &   !< constrained fraction of lateral (horizontal) water loss of 'w_soil_sl_old' (prev. timestep) [unitless]
                                    gw_runoff_sl,     & !< lateral (horizontal) soilwater to groundwater runoff per layer           [kg m-2 yr-1]
                                    percolation_sl,   & !< fraction of material transferred to below layer                          [fraction s-1]
                                    drainage_sl,      & !< vertical water flow (drainage) per soil layer across soil layers         [kg m-2 s-1]
                                    w_soil_sl,        & !< Water content in soil layers                                             [m]
                                    w_soil_sat_sl,    & !< Water content of soil layers at saturation
                                                        !! (derived from volumetric capacity at saturation)                         [m]
                                    w_soil_fc_sl,     & !< Water content of soil layers at field capacity
                                                        !! (derived from volumetric field capacity)                                 [m]
                                    w_soil_pwp_sl,    & !< Water content of soil layers at permanent wilting point
                                                        !! (derived from vol. perm wilt point)                                      [m]
                                    saxtonA,          & !< coefficient in moisture-tension relationship [unitless]
                                    saxtonB,          & !< coefficient in moisture-tension relationship [unitless]
                                    saxtonC,          & !< exponent of moisture-conductivity relationship [unitless]
                                    kdiff_sat_sl,     & !< saturated hydraulic conductivity                                         [m s-1]
                                    w_ice_sl            !< ice content of soil layer                                                [m]

    ! heat fluxes
    TYPE(t_jsb_var_real2d)          :: &
                                    sensible_heat_flx     , &  !< sensible heat flux [W m-2]
                                    latent_heat_flx       , &  !< latent heat flux [W m-2]
                                    ground_heat_flx       , &  !< ground heat flux [W m-2]
                                    ground_heat_flx_old        !< ground heat flux of previous timestep [W m-2]

    ! transpiration
    TYPE(t_jsb_var_real2d)          :: &
                                    transpiration       !< transpiration in [kg m-2 s-1]

    !Soil Water addons for SLM in HYDRO
    TYPE(t_jsb_var_real2d)          :: &
                                    w_soil_root_theta, & !< root zone integrated plant available water [fraction of maximum saturation]
                                    w_soil_root_pot      !< soil water potential [MPa]

    TYPE(t_jsb_var_real3d)          :: &
                                    w_soil_pot_sl        !< soil water potential per layer [MPa]

    ! soil energy
    TYPE(t_jsb_var_real3d)          :: &
                                    t_soil_sl            !< soil temperature [K]

    TYPE(t_jsb_var_real3d)          :: &
      & w_soil_freeze_flux, &                            !< soil water flux from liquid water to ice       [m]
      & w_soil_melt_flux                                 !< soil water flux from ice to liquid water       [m]

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

    ! variables for snow
    TYPE(t_jsb_var_real2d)          :: snow_height, &           !< height (thickness) of all snow layers                         [m]
                                       snow_soil_heat_flux, &   !< heat flux between snow and soil                          [W m-2]
                                       snow_srf_heat_flux,  &   !< heat flux between atmoshere and snow                     [W m-2]
                                       snow_melt_to_soil        !< Snow melt water to soil flux                              [kg m-2 s-1]

    TYPE(t_jsb_var_real3d)          :: snow_present_snl       , &  !< variable to check if there is any snow for each layer     [unitless]
                                       t_snow_snl             , &  !< temperature of snow layer                                 [K]
                                       w_snow_snl             , &  !< water content of snow layer                               [m]
                                       snow_lay_thickness_snl , &  !< snow layer thickness                                      [m]
                                       w_snow_max_snl         , &  !< maximum water content of snow layer                       [m]
                                       w_liquid_snl                !< Liquid water in snow layers                               [m]

    !> ---------------------------------------------------------------------------------------------------- !
    !> temporary SPQ_ variables - used during implementation of using jsbach physics in quincy
    !>

    !> SEB_
    !>
    ! Common variables
    TYPE(t_jsb_var_real2d) ::   &
      & t,                      & !< Surface temperature                                             [K]
      & t_old,                  & !< Surface temperature at previous timestep                        [K]
      & t_unfilt,               & !< Surface temperature, unfiltered (before Asselin)                [K]
      & t_rad4,                 & !< Surface radiative temperature **4                               [K**4]
      & t_eff4,                 & !< Effective surface temperature **4                               [K**4]
                                  !! This is the surface temperature used in the atmospheric
                                  !! radiative heating
      ! & qsat_star,              & !< Saturation specific humidity at surface (qsat^star, see manual) [kg kg-1]
      !                             !! This is the qsat_srf that is used in the vertical diffusion
      !                             !! scheme in the atmosphere
      & heat_cap,               & !< Heat capacity of surface layer                                  [J m-2 K-1]
      & latent_hflx,            & !< Latent heat flux                                                [W m-2]
      & sensible_hflx,          & !< Sensible heat flux                                              [W m-2]
      & forc_hflx,              & !< Additional heat flux at surface                                 [W m-2]
                                  !! (additional to radiative fluxes and latent and sensible
                                  !!  heat fluxes, e.g. ground heat flux)
      & richardson                !< Richardson number (standalone model)                            []


    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & t_lwtr,               & !< Lake surface temperature (water)                                [K]
      & t_lice,               & !< Lake surface temperature (ice)                                  [K]
      & s_lwtr,               & !< Surface dry static energy over lake water                       [m2 s-2]
      & qsat_lwtr,            & !< Saturation specific humidity at surface over lake water        [kg kg-1]
      & s_lice,               & !< Surface dry static energy over lake ice                         [m2 s-2]
      & qsat_lice,            & !< Saturation specific humidity at surface over lake ice           [kg kg-1]
      & latent_hflx_wtr,      &
      & sensible_hflx_wtr,    &
      & latent_hflx_ice,      &
      & sensible_hflx_ice,    &
      & depth_lice,           & !< Thickness of ice on lake                                        [m]
      & fract_lice,           & !< Fraction of lake ice                                            []
      & hflx_form_ice,        & !< Heat flux for ice formation                                     [W m-2]
      & hflx_cond_ice,        & !< Conductive heat flux for ice on lakes
      & hflx_melt_lice,       & !< Heat flux from melting ice on lakes
      & hflx_melt_snow          !< Heat flux from melting snow on lake ice

    !> HYDRO_
    !>
    TYPE(t_jsb_var_real2d) :: &
      & fract_snow              !< Snow fraction of tile                                                []
                                !! 'fract_snow(ic) = fract_snow_soil(ic)' in HYDRO_ interface

    TYPE(t_jsb_var_real2d) :: &
      & fract_wet               !< Wet (skin and pond reservoir) fraction of tile (soil and canopy)         []

    TYPE(t_jsb_var_real2d) :: &
      & fract_pond,           & !< Fraction of tile assumed inundated at given water level              []
      & wtr_pond,             & !< Liquid water content of pond reservoir                               [m water equivalent]
      & ice_pond                !< Frozen water content of pond reservoir                               [m water equivalent]

    ! Additional variables for land type
    TYPE(t_jsb_var_real2d) :: &
      & fract_snow_can,       & !< Snow fraction on canopy                                             []
      & fract_snow_soil,      & !< Snow fraction on soil                                               []
      & weq_snow_soil,        & !< Water content of snow reservoir of soil                             [m water equivalent]
      & snow_soil_dens,       & !< Density of snow on soil                                             [kg m-3]
      & evapotrans_lnd,       & !< Evapotranspiration over vegetation and bare soil type               [kg m-2 s-1]
      & q_snocpymlt             !< Heating by melting of snow on canopy                                [W m-2]

    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & evapo_wtr,            & !< Evaporation over lake water
      & evapo_ice,            & !< Evaporation over lake ice
      & fract_snow_lice,      & !< Fraction of swow on lake ice                                             []
      & weq_snow_lice           !< Water content of snow on lake ice                                        [m water equivalent]


  CONTAINS
    PROCEDURE :: Init => Init_spq_memory
  END TYPE t_spq_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SPQ_ process
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_spq_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,         ONLY: BASIC , MEDIUM, FULL
    USE mo_jsb_io,              ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_CONSTANT
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,     ONLY: t_jsb_model
    USE mo_quincy_output_class, ONLY: unitless
    USE mo_jsb_physical_constants, ONLY: &
      & Tzero,                           &
      & dens_snow,                       & ! Density of snow (l_dynsnow=.FALSE.)
      & dens_snow_min                      ! Density of fresh snow (l_dynsnow=.TRUE.)
    dsl4jsb_Use_processes SPQ_
    dsl4jsb_Use_config(SPQ_)
    USE mo_jsb_control,            ONLY: jsbach_runs_standalone
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_spq_memory), INTENT(inout), TARGET  :: mem             !< spq memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !< process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !< tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: lib_id          !< id of primary lct in lctlib
    INTEGER,              INTENT(in)            :: model_id        !< model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SPQ_)
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb                ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_snow_spq               ! Vertical grid
    INTEGER                    :: table
    TYPE(t_grib2)              :: grib2_desc
    REAL(wp)                   :: ini_snow_dens
    LOGICAL                    :: l_ponds
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_spq_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model          => Get_model(model_id)
    table          =  tables(1)
    hgrid          => Get_grid(model%grid_id)
    surface        => Get_vgrid('surface')
    vgrid_soil_sb  => Get_vgrid('soil_layer_sb')
    vgrid_snow_spq => Get_vgrid('snow_layer_spq')

    dsl4jsb_Get_config(SPQ_)

    IF (dsl4jsb_Config(SPQ_)%l_dynsnow) THEN !SSE_)%l_dynsnow) THEN
      ini_snow_dens = dens_snow_min
    ELSE
      ini_snow_dens = dens_snow
    END IF
    l_ponds = dsl4jsb_Config(SPQ_)%l_ponds ! HYDRO_

    ! ----------------------------------------------------------------------------------------------------- !
    ! add memory only for LAND & PFT tiles
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      CALL mem%Add_var('soil_depth', mem%soil_depth, &
        & hgrid, surface, &
        & t_cf('soil_depth', 'm', 'Soil depth derived from textures (bedrock), for hydrology'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('root_depth', mem%root_depth, &
        & hgrid, surface, &
        & t_cf('root_depth', 'm', 'actual rooting depth'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_srf_new', mem%t_srf_new, &
        & hgrid, surface, &
        & t_cf('t_srf_new', 'K', 'current surface temperature'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 273.15_wp)

      CALL mem%Add_var('t_srf_old', mem%t_srf_old, &
        & hgrid, surface, &
        & t_cf('t_srf_old', 'K', 'previous surface temperature'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 273.15_wp)

      CALL mem%Add_var('temp_srf_eff_4', mem%temp_srf_eff_4, &
        & hgrid, surface, &
        & t_cf('temp_srf_eff_4', 'K**4.0', 'effective surface temperature ** 4.0'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 273.15_wp ** 4.0_wp)

      CALL mem%Add_var('zril_old', mem%zril_old, &
        & hgrid, surface, &
        & t_cf('zril_old', '', 'previous timesteps Reynolds number'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('elevation', mem%elevation, &
        & hgrid, surface, &
        & t_cf('elevation', 'm', 'Mean orography'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart=.FALSE., &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT)

      CALL mem%Add_var('qsat_star', mem%qsat_star, &
        & hgrid, surface, &
        & t_cf('qsat_star', 'kg kg-1 ?', 'saturation specific humidity at surface'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart=.TRUE., &
        & initval_r=0.0075_wp)

      CALL mem%Add_var('s_star', mem%s_star, &
        & hgrid, surface, &
        & t_cf('s_star', 'm2 s-2 ?', 'surface dry static energy'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart=.TRUE., &
        & initval_r=2.9E5_wp)

      CALL mem%Add_var('evapotranspiration', mem%evapotranspiration, &
        & hgrid, surface, &
        & t_cf('evapotranspiration', 'kg m-2 s-1 ?', 'evapotranspiration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart=.TRUE., &
        & initval_r=0.0_wp)

      CALL mem%Add_var('num_sl_above_bedrock', mem%num_sl_above_bedrock, &
        & hgrid, surface, &
        & t_cf('num_sl_above_bedrock', '[unitless]', 'number of soil layers above bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart=.TRUE., &
        & initval_r=0.0_wp)

      CALL mem%Add_var('q_sat_srf', mem%q_sat_srf, &
        & hgrid, surface, &
        & t_cf('q_sat_srf', '[?]', 'surface saturated specific humidity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart=.TRUE., &
        & initval_r=0.0_wp)

      CALL mem%Add_var('soil_depth_sl', mem%soil_depth_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_depth_sl', 'm', 'Width of each soil layer that can be saturated with water, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_lay_width_sl', mem%soil_lay_width_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_lay_width_sl', 'm', 'Width of each soil layer that can be saturated with water, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_lay_depth_center_sl', mem%soil_lay_depth_center_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_lay_depth_center_sl', 'm', 'Depth at the center of each soil layer, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_lay_depth_ubound_sl', mem%soil_lay_depth_ubound_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_lay_depth_ubound_sl', 'm', 'Depth at the upper boundary of each soil layer, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('soil_lay_depth_lbound_sl', mem%soil_lay_depth_lbound_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('soil_lay_depth_lbound_sl', 'm', 'Depth at the lower boundary of each soil layer, until bedrock'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('interception', mem%interception, &
        & hgrid, surface, &
        & t_cf('interception', 'kg m-2 s-1', 'surface interception evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('evapopot', mem%evapopot, &
        & hgrid, surface, &
        & t_cf('evapopot', 'kg m-2 s-1', 'potential evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('evaporation', mem%evaporation, &
        & hgrid, surface, &
        & t_cf('evaporation', 'kg m-2 s-1', 'surface evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('evaporation_snow', mem%evaporation_snow, &
        & hgrid, surface, &
        & t_cf('evaporation_snow', 'kg m-2 s-1', 'snow surface evaporation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('srf_runoff', mem%srf_runoff, &
        & hgrid, surface, &
        & t_cf('srf_runoff', 'kg m-2 s-1', 'surface runoff'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('drainage', mem%drainage, &
        & hgrid, surface, &
        & t_cf('drainage', 'kg m-2 s-1', 'drainage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('gw_runoff', mem%gw_runoff, &
        & hgrid, surface, &
        & t_cf('gw_runoff', 'kg m-2 s-1', 'sum of lateral runoff across all soil layers'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('drainage_fraction', mem%drainage_fraction, &
        & hgrid, surface, &
        & t_cf('drainage_fraction', '1 s-1', 'water mass fraction lost from the soil column by drainage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_skin', mem%w_skin, &
        & hgrid, surface, &
        & t_cf('w_skin', 'm', 'Water content in skin reservoir'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root', mem%w_soil_root, &
        & hgrid, surface, &
        & t_cf('w_soil_root', 'm', 'Water content in root zone of the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root_fc', mem%w_soil_root_fc, &
        & hgrid, surface, &
        & t_cf('w_soil_root_fc', 'm', 'Water content at field capacity in root zone of the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root_pwp', mem%w_soil_root_pwp, &
        & hgrid, surface, &
        & t_cf('w_soil_root_pwp', 'm', 'Water content at permanent wilting point in root zone of the soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('frac_w_lat_loss_sl', mem%frac_w_lat_loss_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('frac_w_lat_loss_sl', '', 'fraction of lateral (horizontal) water loss of w_soil_sl_old'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('gw_runoff_sl', mem%gw_runoff_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('gw_runoff_sl', 'kg m-2 s-1', 'lateral groundwater runoff per layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('percolation_sl', mem%percolation_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('percolation_sl', 'fraction s-1', 'fraction of mass transferred to below layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('drainage_sl', mem%drainage_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('drainage_sl', 'kg m-2 s-1', 'drainage of soil layers'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_sl', mem%w_soil_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_sl', 'm', 'Water content in soil layers'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_sat_sl', mem%w_soil_sat_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_sat_sl', 'm', 'Water content of soil layers at saturation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_fc_sl', mem%w_soil_fc_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_fc_sl', 'm', 'Water content of soil layers at field capacity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_pwp_sl', mem%w_soil_pwp_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_pwp_sl', 'm', 'Water content of soil layers at permanent wilting point'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('saxtonA', mem%saxtonA, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('saxtonA', unitless, 'coefficient in moisture-tension relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('saxtonB', mem%saxtonB, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('saxtonB', unitless, 'coefficient in moisture-tension relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('saxtonC', mem%saxtonC, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('saxtonC', unitless, 'exponent of moisture-conductivity relationship'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('kdiff_sat_sl', mem%kdiff_sat_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('kdiff_sat_sl', 'm s-1', 'saturated hydraulic conductivity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_ice_sl', mem%w_ice_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_ice_sl', 'm', 'Soil layer ice content'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sensible_heat_flx', mem%sensible_heat_flx, &
        & hgrid, surface, &
        & t_cf('sensible_heat_flx', 'W m-2', 'sensible heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('latent_heat_flx', mem%latent_heat_flx, &
        & hgrid, surface, &
        & t_cf('latent_heat_flx', 'W m-2', 'latent heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ground_heat_flx', mem%ground_heat_flx, &
        & hgrid, surface, &
        & t_cf('ground_heat_flx', 'W m-2', 'ground heat flux'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('ground_heat_flx_old', mem%ground_heat_flx_old, &
        & hgrid, surface, &
        & t_cf('ground_heat_flx_old', 'W m-2', 'ground heat flux of previous timestep'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('transpiration', mem%transpiration, &
        & hgrid, surface, &
        & t_cf('transpiration', 'kg m-2 s-1', 'transpiration in'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root_theta', mem%w_soil_root_theta, &
        & hgrid, surface, &
        & t_cf('w_soil_root_theta', 'fraction of maximum saturation', 'root zone integrated plant avail water, fract of max'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_root_pot', mem%w_soil_root_pot, &
        & hgrid, surface, &
        & t_cf('w_soil_root_pot', 'MPa', 'soil water potential'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_pot_sl', mem%w_soil_pot_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_pot_sl', 'MPa', 'soil water potential per layer'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('t_soil_sl', mem%t_soil_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('t_soil_sl', 'K', 'soil temperature'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 273.15_wp)

      CALL mem%Add_var('w_soil_freeze_flux', mem%w_soil_freeze_flux, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_freeze_flux', '-', 'Water flux from liquid to ice from freezing'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('w_soil_melt_flux', mem%w_soil_melt_flux, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('w_soil_melt_flux', '-', 'Water transfer from ice to liquid from melting'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('heat_capa_sl', mem%heat_capa_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('heat_capa_sl', 'J kg-1', 'soil heat capacity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('bulk_dens_sl', mem%bulk_dens_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('bulk_dens_sl', 'kg m-3', 'soil bulk density'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('therm_cond_sl', mem%therm_cond_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('therm_cond_sl', 'W m-2 K-1', 'soil thermal conductivity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('sand_sl', mem%sand_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('sand_sl', '', 'soil sand proportion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('silt_sl', mem%silt_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('silt_sl', '', 'soil silt proportion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('clay_sl', mem%clay_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('clay_sl', '', 'soil clay proportion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('volume_min_sl', mem%volume_min_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('volume_min_sl', '', 'soil mineral volumetric proportion'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_org_min_sl', mem%qmax_org_min_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_org_min_sl', 'mol C kg-1 mineral soil', 'maximum OM sorption capacity of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_po4_min_sl', mem%qmax_po4_min_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_po4_min_sl', 'mol P kg-1 mineral soil', 'maximum po4 sorption capacity of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_nh4_min_sl', mem%qmax_nh4_min_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_nh4_min_sl', 'mol N kg-1 mineral soil', 'maximum po4 sorption capacity of mineral soil'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('qmax_po4_om_sl', mem%qmax_po4_om_sl, &
        & hgrid, vgrid_soil_sb, &
        & t_cf('qmax_po4_om_sl', 'mol P kg-1 SOM', 'maximum po4 sorption capacity of SOM'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('snow_height', mem%snow_height, &
          &  hgrid, surface, &
          &  t_cf('snow_height','m', 'snow depth'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .TRUE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var( 'snow_soil_heat_flux', mem%snow_soil_heat_flux, &
          &  hgrid, surface, &
          &  t_cf('snow_soil_heat_flux','W m-2', 'snow soil heat flux'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('snow_srf_heat_flux', mem%snow_srf_heat_flux, &
          &  hgrid, surface, &
          &  t_cf('snow_srf_heat_flux','W m-2', 'snow srf heat flux'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('snow_melt_to_soil', mem%snow_melt_to_soil, &
          &  hgrid, surface, &
          &  t_cf('snow_melt_to_soil','kg m-2 s-1', 'snow melt water to soil'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('snow_present_snl', mem%snow_present_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('snow_present_snl',unitless, 'check if there is snow'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('t_snow_snl', mem%t_snow_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('t_snow_snl','m', 'temperature snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('w_snow_snl', mem%w_snow_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_snow_snl','m', 'water content snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('snow_lay_thickness_snl', mem%snow_lay_thickness_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('snow_lay_thickness_snl','m', 'snow layer thickness'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('w_snow_max_snl', mem%w_snow_max_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_snow_max_snl','m', 'maximum water content snow layer'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

      CALL mem%Add_var('w_liquid_snl', mem%w_liquid_snl, &
          &  hgrid, vgrid_snow_spq, &
          &  t_cf('w_liquid_snl','m', 'w_liquid_snl'), &
          &  t_grib1(table, 255, grib_bits), &
          &  t_grib2(255, 255, 255, grib_bits), &
          &  prefix, suffix, &
          & output_level = FULL, &
          &  loutput = .FALSE., &
          &  lrestart = .TRUE., &
          &  initval_r = 0.0_wp) ! initvalue

    END IF ! IF One_of(LAND_TYPE  OR VEG_TYPE, lct_ids)

    !> ---------------------------------------------------------------------------------------------------- !
    !> temporary SPQ_ variables - used during implementation of using jsbach physics in quincy
    !>

    CALL mem%Add_var( 't', mem%t,                                                    &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp', 'K', 'surface temperature'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC, &
      & initval_r=280.0_wp )

    ! IF (model%config%use_tmx) THEN
    !   CALL mem%Add_var( 't_rad**4', mem%t_rad4,                                      &
    !     & hgrid, surface,                                                            &
    !     & t_cf('sfc_temp_rad**4', 'K', 'surface radiative temperature **4'),         &
    !     & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
    !     & prefix, suffix,                                                            &
    !     & output_level=BASIC,                                                        &
    !     & initval_r=280.0_wp**4 )
    ! END IF

    CALL mem%Add_var( 't_old', mem%t_old,                                            &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_old', 'K', 'surface temperature at previous timestep'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & initval_r=280.0_wp )

    CALL mem%Add_var( 't_unfilt', mem%t_unfilt,                                      &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_unfilt', 'K', 'unfiltered surface temperature'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=280.0_wp )

    CALL mem%Add_var( 'heat_cap', mem%heat_cap,                                      &
      & hgrid, surface,                                                              &
      & t_cf('sfc_heat_cap', 'J m-2 K-1)', 'surface layer heat capacity'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & output_level=BASIC, &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'latent_hflx', mem%latent_hflx,                                &
      & hgrid, surface,                                                              &
      & t_cf('latent_hflx', 'W m-2', 'latent heat flux density at surface'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'sensible_hflx', mem%sensible_hflx,                            &
      & hgrid, surface,                                                              &
      & t_cf('sensible_hflx', 'W m-2', 'sensible heat flux at surface'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'forc_hflx', mem%forc_hflx,                                    &
      & hgrid, surface,                                                              &
      & t_cf('forc_hflx', 'W m-2)', 'additional heat flux at surface'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & lrestart=.FALSE.,                                                            &
      & output_level=BASIC, &
      & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    IF (jsbach_runs_standalone()) THEN
      CALL mem%Add_var( 'richardson', mem%richardson,                                &
        & hgrid, surface,                                                            &
        & t_cf('richardson', '-', 'Richardson number'),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.TRUE.,                                                           &
        & initval_r=0.0_wp )
    END IF

    ! IF (model%config%use_tmx .OR. One_of(LAND_TYPE, lct_ids(:)) > 0) THEN
    !   CALL mem%Add_var( 't_eff**4', mem%t_eff4,                                              &
    !     & hgrid, surface,                                                                    &
    !     & t_cf('surface_temperature_eff**4', 'K', 'effective surface temperature**4'),       &
    !     & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
    !     & prefix, suffix,                                                                    &
    !     & lrestart=.FALSE.,                                                                  &
    !     & initval_r=280.0_wp**4 )
    ! END IF

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 't_lwtr', mem%t_lwtr,                                          &
        & hgrid, surface,                                                              &
        & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr '),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=280.0_wp )

      CALL mem%Add_var( 's_lwtr', mem%s_lwtr,                                                  &
        & hgrid, surface,                                                                      &
        & t_cf('sfc_dry_static_energy_wtr', 'm2 s-2', 'surface dry static energy (lake wtr)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & lrestart=.FALSE.,                                                                    &
        & initval_r=2.9E5_wp )

      CALL mem%Add_var( 'qsat_lwtr', mem%qsat_lwtr,                                              &
        & hgrid, surface,                                                                        &
        & t_cf('sfc_qsat_lwtr', 'm2 s-2', 'surface specific humidity at saturation (lake wtr)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & lrestart=.FALSE.,                                                                      &
        & initval_r=0.0075_wp )

      CALL mem%Add_var( 'latent_hflx_wtr', mem%latent_hflx_wtr,                        &
        & hgrid, surface,                                                              &
        & t_cf('latent_hflx_wtr', 'W m-2', 'latent heat flux over lake water'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE.,                                                             &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'sensible_hflx_wtr', mem%sensible_hflx_wtr,                    &
        & hgrid, surface,                                                              &
        & t_cf('sensible_hflx_wtr', 'W m-2', 'sensible heat flux over lake water'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE.,                                                             &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! Fraction of ice on lake
      CALL mem%Add_var( 'fract_lice', mem%fract_lice,                                  &
        & hgrid, surface,                                                              &
        & t_cf('fract_lice', '', 'fraction of ice on lake'),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      IF (dsl4jsb_Config(SPQ_)%l_ice_on_lakes) THEN !SEB_)%l_ice_on_lakes) THEN

        CALL mem%Add_var( 't_lice', mem%t_lice,                                          &
          & hgrid, surface,                                                              &
          & t_cf('sfc_temp_ice', 'K', 'lake surface temperature ice '),                  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=Tzero )

        CALL mem%Add_var( 's_lice', mem%s_lice,                                                  &
          & hgrid, surface,                                                                      &
          & t_cf('sfc_dry_static_energy_ice', 'm2 s-2', 'surface dry static energy (lake ice)'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
          & prefix, suffix,                                                                      &
          & lrestart=.FALSE.,                                                                    &
          & initval_r=2.9E5_wp )

        CALL mem%Add_var( 'qsat_lice', mem%qsat_lice,                                              &
          & hgrid, surface,                                                                        &
          & t_cf('sfc_qsat_lice', 'm2 s-2', 'surface specific humidity at saturation (lake ice)'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix,                                                                        &
          & lrestart=.FALSE.,                                                                      &
          & initval_r=0.0075_wp )

        CALL mem%Add_var( 'latent_hflx_ice', mem%latent_hflx_ice,                        &
          & hgrid, surface,                                                              &
          & t_cf('latent_hflx_ice', 'W m-2', 'latent heat flux over lake ice'),          &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & lrestart=.TRUE.,                                                             &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'sensible_hflx_ice', mem%sensible_hflx_ice,                    &
          & hgrid, surface,                                                              &
          & t_cf('sensible_hflx_ice', 'W m-2', 'sensible heat flux over lake ice'),      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & lrestart=.TRUE.,                                                             &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'depth_lice', mem%depth_lice,                                  &
          & hgrid, surface,                                                              &
          & t_cf('ice_thickness', 'm', 'thickness of ice on lake'),                      &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp )

        CALL mem%Add_var( 'hflx_form_ice', mem%hflx_form_ice,                            &
          & hgrid, surface,                                                              &
          & t_cf('hflx_form_ice', 'W m-2)', 'heat flux for ice formation on lake ice'),  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'hflx_cond_ice', mem%hflx_cond_ice,                            &
          & hgrid, surface,                                                              &
          & t_cf('hflx_cond_ice', 'W m-2)', 'conductive heat flux on lake ice'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'hflx_melt_lice', mem%hflx_melt_lice,                          &
          & hgrid, surface,                                                              &
          & t_cf('hflx_melt_lice', 'W m-2)', 'heat flux from ice melt on lakes'),        &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

        CALL mem%Add_var( 'hflx_melt_snow', mem%hflx_melt_snow,                          &
          & hgrid, surface,                                                              &
          & t_cf('hflx_melt_snow', 'W m-2)', 'heat_flux from snow melt on lakes'),       &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
          & prefix, suffix,                                                              &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      END IF  ! l_ice_on_lakes
    END IF    ! LAKE_TYPE


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

    IF ( (     One_of(VEG_TYPE,     lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE,    lct_ids(:)) > 0 &
      &   .OR. One_of(GLACIER_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE,    lct_ids(:)) > 0 &
      &  ) ) THEN

      IF (TRIM(mem%owner_tile_name) == 'box') THEN
        grib2_desc = t_grib2(1,0,201, grib_bits)
      ELSE
        grib2_desc = t_grib2(255, 255, 255, grib_bits)
      END IF
      CALL mem%Add_var( 'fract_wet', mem%fract_wet,                              &
        & hgrid, surface,                                                        &
        & t_cf('fract_wet', '-', 'Wet surface fraction'),                        &
        & t_grib1(table, 255, grib_bits), grib2_desc,                            &
        & prefix, suffix,                                                        &
        & output_level=BASIC,                                                    &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'fract_pond', mem%fract_pond,                                  &
        & hgrid, surface,                                                              &
        & t_cf('fract_pond', '-', 'Inundated surface fraction'),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=l_ponds, lrestart_cont=.TRUE.       ,                               &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'wtr_pond', mem%wtr_pond,                                      &
        & hgrid, surface,                                                              &
        & t_cf('wtr_pond', 'm (water equivalent)', 'Liquid water in pond reservoir'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=l_ponds, lrestart_cont=.TRUE.       ,                               &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'ice_pond', mem%ice_pond,                                      &
        & hgrid, surface,                                                              &
        & t_cf('ice_pond', 'm (water equivalent)', 'Frozen water in pond reservoir'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=l_ponds, lrestart_cont=.TRUE.       ,                               &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'fract_snow_can', mem%fract_snow_can,                        &
        & hgrid, surface,                                                            &
        & t_cf('fract_snow_can', '-', 'Snow fraction on canopy'),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'fract_snow_soil', mem%fract_snow_soil,                      &
        & hgrid, surface,                                                            &
        & t_cf('fract_snow_soil', '-', 'Snow fraction on soil'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL mem%Add_var( 'weq_snow_soil', mem%weq_snow_soil,                          &
        & hgrid, surface,                                                            &
        & t_cf('weq_snow_soil', 'm (water equivalent)', 'Snow amount on soil'),      &
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
      CALL mem%Add_var( 'evapotrans_lnd', mem%evapotrans_lnd,                                                     &
        & hgrid, surface,                                                                                         &
        & t_cf('surface_evapotranspiration_lnd', 'kg m-2 s-1', 'Evapotranspiration from land surface w/o lakes'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & lrestart=.FALSE.,                                                                                       &
        & loutput=.TRUE., output_level=BASIC,                                                                     &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      ! IF (.NOT. model%config%use_tmx) THEN
      !   CALL mem%Add_var( 'q_snocpymlt', mem%q_snocpymlt,                                 &
      !     & hgrid, surface,                                                               &
      !     & t_cf('heating_snow_cpy_melt', 'W m-2', 'Heating due to snow melt on canopy'), &
      !     & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
      !     & prefix, suffix,                                                               &
      !     & lrestart=.FALSE.,                                                             &
      !     & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      ! END IF
    END IF   ! lct

    ! Additional variables for LAKE lct
    IF (      One_of(LAKE_TYPE, lct_ids(:)) > 0 ) THEN

      CALL mem%Add_var( 'evapo_wtr', mem%evapo_wtr,                                &
        & hgrid, surface,                                                          &
        & t_cf('evapo_wtr', 'kg m-2 s-1', 'Evaporation from lake water'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                          &
        & lrestart=.FALSE.,                                                        &
        & initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (dsl4jsb_Config(SPQ_)%l_ice_on_lakes) THEN !SEB_)%l_ice_on_lakes) THEN
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

        CALL mem%Add_var( 'weq_snow_lice', mem%weq_snow_lice,                        &
          & hgrid, surface,                                                          &
          & t_cf('weq_snow_lice', 'm water equivalent', 'Snow amount on lake ice'),  &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
          & prefix, suffix,                                                          &
          & output_level=BASIC,                                                      &
          & initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      END IF
    END IF

  END SUBROUTINE Init_spq_memory

#endif
END MODULE mo_spq_memory_class
