!>
!! @brief Contains methods for the atmosphere-land interface memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-27)
!!
MODULE mo_atmland_memory

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory, t_land2atm_memory
  USE mo_cdi_constants,     ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, &
       ZA_SURFACE, &
       FILETYPE_NC2, &
       GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_atm2land_memory, new_land2atm_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_memory'

CONTAINS

  FUNCTION new_atm2land_memory(varlist_name, prefix, suffix, tile) RESULT(atm2land)

    USE mo_jsb_varlist_iface, ONLY: t_var_list, get_var_list, new_var_list, default_var_list_settings, add_var
    USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_base,          ONLY: get_model

    CHARACTER(len=*),           INTENT(in) :: varlist_name
    CHARACTER(len=*),           INTENT(in) :: prefix
    CHARACTER(len=*),           INTENT(in) :: suffix
#ifdef HAVE_F2003
    TYPE(t_jsb_tile), POINTER , INTENT(in) :: tile
#else
    TYPE(t_jsb_tile), POINTER              :: tile
#endif
    TYPE(t_atm2land_memory), POINTER       :: atm2land

    TYPE(t_var_list),        POINTER :: this_list
    TYPE(t_jsb_model),       POINTER :: model

    INTEGER :: nproma, nblks
    INTEGER :: shape2d(2)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_atm2land_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(atm2land)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL new_var_list(this_list, TRIM(varlist_name), patch_id=model%grid%patch%id)
      CALL default_var_list_settings(this_list, lrestart=.TRUE., restart_type=FILETYPE_NC2)
    END IF

    nproma       = model%grid%nproma
    nblks        = model%grid%nblks

    shape2d   = (/ nproma, nblks /)

    CALL add_var( this_list, prefix//'t_air'//suffix, atm2land%t_air,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('air_temperature', 'K', '', DATATYPE_FLT32),                           &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'q_air'//suffix, atm2land%q_air,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('air_specific_humidity', '', '', DATATYPE_FLT32),                      &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'rain'//suffix, atm2land%rain,                        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('rain_fall', 'kg/m2/s', '', DATATYPE_FLT32),                           &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'snow'//suffix, atm2land%snow,                        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('snow_fall', 'kg/m2/s', '', DATATYPE_FLT32),                           &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'wind_air'//suffix, atm2land%wind_air,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('air_wind_speed', 'm/s', '', DATATYPE_FLT32),                          &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'wind_10m'//suffix, atm2land%wind_10m,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('10m_wind_speed', 'm/s', '', DATATYPE_FLT32),                          &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'lwrad_srf_down'//suffix, atm2land%lwrad_srf_down,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('downward_longwave_radiation', '', '', DATATYPE_FLT32),                &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'swrad_srf_down'//suffix, atm2land%swrad_srf_down,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('downward_shortwave_radiation', 'K', '', DATATYPE_FLT32),              &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'press_srf'//suffix, atm2land%press_srf,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('surface_pressure', '', '', DATATYPE_FLT32),                           &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'drag_srf'//suffix, atm2land%drag_srf,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('surface_drag', '', '', DATATYPE_FLT32),                               &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'t_acoef'//suffix, atm2land%t_acoef,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('temperature_acoef', '', '', DATATYPE_FLT32),                          &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'t_bcoef'//suffix, atm2land%t_bcoef,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('temperature_acoef', '', '', DATATYPE_FLT32),                          &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'q_aceof'//suffix, atm2land%q_acoef,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('specific_humidity_acoef', '', '', DATATYPE_FLT32),                    &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'q_bcoef'//suffix, atm2land%q_bcoef,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('specific_humidity_bcoef', '', '', DATATYPE_FLT32),                    &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'cos_zenith_angle'//suffix, atm2land%cos_zenith_angle,&
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                             &
         & t_cf_var('cosine_zenith_angle', '', '', DATATYPE_FLT32),                        &
         & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),         &
         & ldims=shape2d, initval_r=0.0_wp )

  END FUNCTION new_atm2land_memory

  FUNCTION new_land2atm_memory(varlist_name, prefix, suffix, tile) RESULT(land2atm)

    USE mo_jsb_varlist_iface, ONLY: t_var_list, get_var_list, new_var_list, default_var_list_settings, add_var
    USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_base,          ONLY: get_model

    CHARACTER(len=*),           INTENT(in) :: varlist_name
    CHARACTER(len=*),           INTENT(in) :: prefix
    CHARACTER(len=*),           INTENT(in) :: suffix
#ifdef HAVE_F2003
    TYPE(t_jsb_tile), POINTER , INTENT(in) :: tile
#else
    TYPE(t_jsb_tile), POINTER              :: tile
#endif
    TYPE(t_land2atm_memory), POINTER       :: land2atm

    TYPE(t_var_list),        POINTER :: this_list
    TYPE(t_jsb_model),       POINTER :: model

    INTEGER :: nproma, nblks
    INTEGER :: shape2d(2)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_land2atm_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(land2atm)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL new_var_list(this_list, TRIM(varlist_name), patch_id=model%grid%patch%id)
      CALL default_var_list_settings(this_list, lrestart=.TRUE., restart_type=FILETYPE_NC2)
    END IF

    nproma       = model%grid%nproma
    nblks        = model%grid%nblks

    shape2d   = (/ nproma, nblks /)

    CALL add_var( this_list, prefix//'t_srf'//suffix, land2atm%t_srf,                      &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('surface_temperature', 'K', '', DATATYPE_FLT32),                          &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'t_rad_srf'//suffix, land2atm%t_rad_srf,              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('radiative_surface_temperature', 'K', '', DATATYPE_FLT32),                &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'qsat_srf'//suffix, land2atm%qsat_srf,                &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('saturated_surface_specific_humidity', '', '', DATATYPE_FLT32),           &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'fact_q_air'//suffix, land2atm%fact_q_air,            &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('factor_air_moisture', '', '', DATATYPE_FLT32),                           &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'fact_qsat_srf'//suffix, land2atm%fact_qsat_srf,      &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('factor_qsat_srf', '', '', DATATYPE_FLT32),                               &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'zh_srf'//suffix, land2atm%zh_srf,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('surface_roughness_humidity', '', '', DATATYPE_FLT32),                    &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'zm_srf'//suffix, land2atm%zm_srf,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('surface_roughness_momentum', '', '', DATATYPE_FLT32),                    &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'alb_vis_dif'//suffix, land2atm%alb_vis_dif,          &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('albedo_visible_diffuse', '', '', DATATYPE_FLT32),                        &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'alb_nir_dif'//suffix, land2atm%alb_nir_dif,          &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('albedo_near_infrared_diffuse', '', '', DATATYPE_FLT32),                  &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'alb_vis_dir'//suffix, land2atm%alb_vis_dir,          &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('albedo_visible_direct', '', '', DATATYPE_FLT32),                         &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

    CALL add_var( this_list, prefix//'alb_nir_dir'//suffix, land2atm%alb_nir_dir,          &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                &
      & t_cf_var('albedo_near_infrared_dir', '', '', DATATYPE_FLT32),                      &
      & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),            &
      & ldims=shape2d, initval_r=0.0_wp )

  END FUNCTION new_land2atm_memory

END MODULE mo_atmland_memory
