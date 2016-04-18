!>
!! @brief Contains the memory class for hydrology.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-29)
!!

MODULE mo_hydro_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,    ONLY: LAND_TYPE, BARE_TYPE, VEG_TYPE, LAKE_TYPE, GLACIER_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hydro_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_hydro_memory

    TYPE(t_jsb_var_real2d) :: &
      & fract_water_srf,  & !< Wet (skin reservoir) fraction of tile                                    []
      & fract_snow_srf,   & !< Snow fraction of tile                                                    []
      & w_skin_srf,       & !< Water content in skin reservoir                                          [m]
      & w_snow_srf,       & !< Water content of snow reservoir                                          [m water equivalent]
      & tte_corr            !< Correction term of TTE for snow melt

    TYPE(t_jsb_var_real2d) :: &
      & evapopot,         & !< Potential evaporation
      & evapotrans          !< Evapotranspiration

    TYPE(t_jsb_var_real2d) :: &
      & runoff,           & !< surface runoff
      & drainage,         & !< drainage
      & discharge,        & !< discharge (local)
      & discharge_ocean,  & !< discharge to the ocean
      & internal_drain      !< internal drain

    ! Additional variables for soil
    TYPE(t_jsb_var_real2d) :: &
      & fract_water_soil, & !< Wet (skin reservoir) fraction of soil                                    []
      & fract_snow_soil,  & !< Snow fraction of soil                                                    []
      & w_skin_soil,      & !< Water content in skin reservoir on soil                                  [m]
      & w_snow_soil,      & !< Water content of snow reservoir of soil                                  [m water equivalent]
      & w_soil_overflow,  & !< Water content of reservoir for w_soil_overflown soil water (terraplanet) [m]
      & w_soil_column       !< Water content in the whole soil column                                   [m]

    TYPE(t_jsb_var_real3d) :: &
      & w_soil_sl,        & !< Water content in soil layers                                             [m]
      & w_soil_sat_sl,    & !< Water content of soil layers at saturation (derived from soil porosity)  [m]
      & w_soil_fc_sl,     & !< Water content of soil layers at field capacity
                            !! (derived from volumetric field capacity)                                 [m]
      & w_soil_pwp_sl       !< Water content of soil layers at permanent wilting point
                            !! (derived from vol. perm wilt point)                                      [m]

    ! Additional variables for PFT lct_type
    TYPE(t_jsb_var_real2d) :: &
      & fract_water_can,  & !< Wet (skin reservoir) fraction of canopy                                  []
      & fract_snow_can,   & !< Snow fraction of canopy                                                  []
      & w_snow_can,       & !< Water content of snow reservoir onn canopy                               [m water equivalent]
      & w_skin_can,       & !< Water content in skin reservoir on canopy                                [m]
      & w_soil_root,      & !< Water content in root zone of the soil                                   [m]
      & w_soil_root_fc,   & !< Water content at field capacity in root zone of the soil                 [m]
      & w_soil_root_pwp,  & !< Water content at permanent wilting point in root zone of the soil        [m]
      & water_stress,     & !< Water stress factor of canopy                        (1: no water stress, 0: infinite stress)
      & trans               !< Transpiration

    ! Additional variables for GLACIER lct_type
    TYPE(t_jsb_var_real2d) :: &
      & w_snow_glacier,   & !< Water content of snow on glacier                                         [m water equivalent]
      & pme_glacier         !< P-E on glaciers

    ! Additional variables if no separate glacier lct_type is used and glaciers are treated as part of SOIL lct_type
    TYPE(t_jsb_var_real2d) :: &
      & fract_glac          !< Glacier fraction of tile                                                 []

    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & ice_depth,        & !< Thickness of ice on lake                                                 [m]
      & ice_fract,        & !< Fraction of ice on lake                                                  []
      & fract_snow_ice,   & !< Snow fraction on lake ice                                                []
      & w_snow_ice          !< Water content of snow on (lake) ice                                      [m water equivalent]

  CONTAINS
    PROCEDURE :: Init => Init_hydro_memory
  END TYPE t_hydro_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_memory_class'

CONTAINS

  SUBROUTINE Init_hydro_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_hydro_memory), INTENT(inout) :: mem
    CHARACTER(len=*),      INTENT(in)    :: prefix
    CHARACTER(len=*),      INTENT(in)    :: suffix
    INTEGER,               INTENT(in)    :: lct_ids(:)
    INTEGER,               INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                         ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer, soil_w ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_hydro_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')
    !canopy_layer => Get_vgrid('canopy_layer')
    soil_w       => Get_vgrid('soil_depth_water')

    IF (ASSOCIATED(mem%fract_water_srf%ptr)) RETURN

    CALL jsb_add_var( mem, 'fract_water_srf', mem%fract_water_srf, &
      & hgrid, surface,                                                            &
      & t_cf('fraction of water on surface surface_wet_fraction', '-', ''),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'fract_snow_srf', mem%fract_snow_srf, &
      & hgrid, surface,                                                          &
      & t_cf('fraction of snow on surface', '-', ''),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'w_skin_srf', mem%w_skin_srf,      &
      & hgrid, surface,                                                       &
      & t_cf('Water content in skin reservoir of surface. surface_skin_reservoir', 'm', ''), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'w_snow_srf', mem%w_snow_srf,      &
      & hgrid, surface,                                                       &
      & t_cf('Water content of snow reservoir on surface', 'm', ''), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                              &
      & initval_r=-1.0_wp )    ! -1 used to determine start of experiment

    CALL jsb_add_var( mem, 'tte_corr', mem%tte_corr,                                     &
      & hgrid, surface,                                                                  &
      & t_cf('tte_correction_snowmelt', '', ''),                                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'evapopot', mem%evapopot,                                     &
      & hgrid, surface,                                                                  &
      & t_cf('surface_evaporation_potential', '', 'Potential evaporation from surface'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'evapotrans', mem%evapotrans,                                 &
      & hgrid, surface,                                                                  &
      & t_cf('surface_evapotranspiration', '', 'Evapotranspiration from surface'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'runoff', mem%runoff,                   &
      & hgrid, surface,                                                                    &
      & t_cf('surface_runoff', 'kg m-2 s-1', 'Surface runoff'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'drainage', mem%drainage,               &
      & hgrid, surface,                                                                    &
      & t_cf('drainage', 'kg m-2 s-1', 'drainage'),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'discharge', mem%discharge,             &
      & hgrid, surface,                                                                    &
      & t_cf('discharge', 'kg m-2 s-1', 'local discharge'),                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'discharge_ocean', mem%discharge_ocean, &
      & hgrid, surface,                                                                    &
      & t_cf('discharge_ocean', 'kg m-2 s-1', 'discharge to the ocean'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'internal_drain', mem%internal_drain,           &
      & hgrid, surface,                                                                    &
      & t_cf('internal_drain', 'kg m-2 s-1', 'internal drainage'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    ! Additional variables for soil
    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(BARE_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  )                                       &
      &  .AND. .NOT. ASSOCIATED(mem%fract_water_soil%ptr) ) THEN

      CALL jsb_add_var( mem, 'fract_water_soil', mem%fract_water_soil, &
        & hgrid, surface,                                                              &
        & t_cf('fraction of water on soil', '-', ''),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, 'fract_snow_soil', mem%fract_snow_soil, &
        & hgrid, surface,                                                            &
        & t_cf('fraction of snow on soil', '-', ''),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_skin_soil', mem%w_skin_soil,    &
        & hgrid, surface,                                                       &
        & t_cf('Water content in skin reservoir on soil', 'm', ''),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_snow_soil', mem%w_snow_soil, &
        & hgrid, surface,                                                    &
        & t_cf('Water content of snow reservoir on soil', 'm', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_soil_overflow', mem%w_soil_overflow, &
        & hgrid, surface,                                                            &
        & t_cf('Water content of reservoir for soil water', 'm', ''),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_soil_column', mem%w_soil_column,     &
        & hgrid, surface,                                                            &
        & t_cf('Water content in the whole soil column', 'm', ''),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_soil_sl', mem%w_soil_sl,     &
        & hgrid, soil_w,                                                     &
        & t_cf('Water content in soil layers', 'm', ''),            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! R: initval =0

      CALL jsb_add_var( mem, 'w_soil_sat_sl', mem%w_soil_sat_sl, &
        & hgrid, soil_w,                                                         &
        & t_cf('Water content in soil layers at saturation', 'm', ''),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'w_soil_fc_sl', mem%w_soil_fc_sl,      &
        & hgrid, soil_w,                                                            &
        & t_cf('Water content in soil layers at field capacity', 'm', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_soil_pwp_sl', mem%w_soil_pwp_sl,             &
        & hgrid, soil_w,                                                                     &
        & t_cf('Water content in soil layers at permanent wilting point', 'm', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! R: initval =0

    END IF

    ! Additional variables for vegetation
    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &  )                                       &
      &  .AND. .NOT. ASSOCIATED(mem%fract_water_can%ptr)) THEN

      CALL jsb_add_var( mem, 'fract_water_can', mem%fract_water_can, &
        & hgrid, surface,                                                            &
        & t_cf('fraction of water on canopy', '-', ''),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'fract_snow_can', mem%fract_snow_can, &
        & hgrid, surface,                                                          &
        & t_cf('fraction of snow on canopy', '-', ''),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_skin_can', mem%w_skin_can,      &
        & hgrid, surface,                                                       &
        & t_cf('Water content in skin reservoir on canopy', 'm', ''),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_snow_can', mem%w_snow_can,     &
        & hgrid, surface,                                                      &
        & t_cf('Water content of snow reservoir on canopy', 'm', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_soil_root', mem%w_soil_root,         &
        & hgrid, surface,                                                            &
        & t_cf('Water content in root zone of the soil', 'm', ''),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_soil_root_fc', mem%w_soil_root_fc,            &
        & hgrid, surface,                                                                     &
        & t_cf('Water content at field capacity in root zone of the soil', 'm', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_soil_root_pwp', mem%w_soil_root_pwp,                   &
        & hgrid, surface,                                                                              &
        & t_cf('Water content at permanent wilting point in root zone of the soil', 'm', ''), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'water_stress', mem%water_stress, &
        & hgrid, surface,                                                      &
        & t_cf('Water stress factor of canopy', '-', ''),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'trans', mem%trans,                                           &
        & hgrid, surface,                                                                  &
        & t_cf('surface_transpiration', '', 'Transpiration from surface'),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

    END IF

    ! Additional variables for GLACIER lct
    IF ( (     One_of(GLACIER_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(LAND_TYPE,    lct_ids(:)) > 0 &
      &  )                                          &
      &  .AND. .NOT. ASSOCIATED(mem%pme_glacier%ptr)) THEN

      CALL jsb_add_var( mem, 'w_snow_glacier', mem%w_snow_glacier,  &
        & hgrid, surface,                                                           &
        & t_cf('Water content of snow reservoir on glacier', 'm', ''),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                              &
        & initval_r=20.0_wp )

      CALL jsb_add_var( mem, 'pme_glacier', mem%pme_glacier,         &
        & hgrid, surface,                                                                    &
        & t_cf('pme_glacier', 'kg m-2 s-1', 'P-E on glaciers'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

    END IF

    ! Additional variables for LAKE lct
    IF (      One_of(LAKE_TYPE, lct_ids(:)) > 0 &
      & .AND. .NOT. ASSOCIATED(mem%fract_snow_ice%ptr)) THEN

      CALL jsb_add_var( mem, 'fract_snow_ice', mem%fract_snow_ice, &
        & hgrid, surface,                                                          &
        & t_cf('fraction of snow on lake ice', '-', ''),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),       &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'w_snow_ice', mem%w_snow_ice,  &
        & hgrid, surface,                                                           &
        & t_cf('Water content of snow reservoir on lake ice', 'm', ''),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                              &
        & initval_r=0.0_wp )

    END IF

    END SUBROUTINE Init_hydro_memory
#endif
END MODULE mo_hydro_memory_class
