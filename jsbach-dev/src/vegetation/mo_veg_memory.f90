!>
!! @brief Contains methods for the vegetation memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_memory
#ifndef __NO_JSBACH__

  USE mo_veg_memory_class,  ONLY: t_veg_memory, t_veg_phys_state, t_veg_water_state
  USE mo_cdi_constants,     ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16,   &
                                  ZA_SURFACE, ZA_DEPTH_BELOW_LAND,                   &
                                  FILETYPE_NC2,                                      &
                                  GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL, &
                                  TSTEP_CONSTANT

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_veg_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory'

CONTAINS

  FUNCTION new_veg_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

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
    TYPE(t_veg_memory), POINTER           :: mem

    TYPE(t_veg_memory), POINTER :: new
    TYPE(t_var_list),   POINTER :: this_list
    TYPE(t_jsb_model), POINTER  :: model

    INTEGER :: nproma, nblks, nsoil_water
    INTEGER :: shape2d(2), shape3d_w(3), shape3d_mon(3)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_veg_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(new)

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
    nsoil_water  = model%config%soil%nsoil_water

    shape2d       = (/ nproma, nblks /)
    shape3d_w     = (/ nproma, nsoil_water, nblks /)
    shape3d_mon   = (/ nproma, nblks, 14 /)

    new%phys_state  => new_phys_state()
    new%water_state => new_water_state()

    mem => new

  CONTAINS

    FUNCTION new_phys_state() RESULT(pstate)

      TYPE(t_veg_phys_state), POINTER :: pstate

      ALLOCATE(pstate)

      CALL add_var( this_list, prefix//'lai'//suffix, pstate%lai,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('leaf_area_index', '', '', DATATYPE_FLT32),                                &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'lai_mon'//suffix, pstate%lai_mon,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('monthly_leaf_area_index', '', '', DATATYPE_FLT32),                        &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape3d_mon, initval_r=0.0_wp, loutput=.FALSE., lrestart=.FALSE. )

      CALL add_var( this_list, prefix//'fract_fpc'//suffix, pstate%fract_fpc,                 &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('foliage_projected_cover', '', '', DATATYPE_FLT32),                        &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'fract_fpc_max'//suffix, pstate%fract_fpc_max,         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('maximum_fpc', '', '', DATATYPE_FLT32),                                    &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'fract_forest'//suffix, pstate%fract_forest,           &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('forest_fraction', '', '', DATATYPE_FLT32),                                &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'canopy_cond'//suffix, pstate%canopy_cond,             &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('canopy_conductance', '', '', DATATYPE_FLT32),                             &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'skin_can_depth'//suffix, pstate%skin_can_depth,       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('canopy_skin_depth', 'm', '', DATATYPE_FLT32),                             &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'root_depth'//suffix, pstate%root_depth,               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('rooting_depth', '', '', DATATYPE_FLT32),                                  &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'root_depth_l'//suffix, pstate%root_depth_l,           &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                        &
        & t_cf_var('root_depth_layer', 'm', '', DATATYPE_FLT32),                              &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape3d_w, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'alb_vis_can'//suffix, pstate%alb_vis_can,             &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('canopy_albedo_visible', '', '', DATATYPE_FLT32),                          &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'alb_nir_can'//suffix, pstate%alb_nir_can,             &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
        & t_cf_var('canopy_albedo_nir', '', '', DATATYPE_FLT32),                              &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_phys_state

    FUNCTION new_water_state() RESULT(wstate)

      TYPE(t_veg_water_state), POINTER :: wstate

      ALLOCATE(wstate)

      CALL add_var( this_list, prefix//'wfract_can'//suffix, wstate%wfract_can,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('wet_fraction_canopy', '', '', DATATYPE_FLT32),                                      &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'wsr_can'//suffix, wstate%wsr_can,                               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('skin_reservoir_canopy', 'm', '', DATATYPE_FLT32),                                   &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'sfract_can'//suffix, wstate%sfract_can,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('snow_fraction_canopy', '', '', DATATYPE_FLT32),                                     &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'wsn_can'//suffix, wstate%wsn_can,                               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('snow_water_equivalent_canopy', 'm', '', DATATYPE_FLT32),                            &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'ws_root'//suffix, wstate%ws_root,                               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('soil_moisture_root', 'm', '', DATATYPE_FLT32),                                      &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'ws_fc_root'//suffix, wstate%ws_fc_root,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('field_capacity_root', 'm', '', DATATYPE_FLT32),                                     &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'ws_pwp_root'//suffix, wstate%ws_pwp_root,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('wilting_root', 'm', '', DATATYPE_FLT32),                                            &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'water_stress'//suffix, wstate%water_stress,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('water_stress_factor', '', '', DATATYPE_FLT32),                                      &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_water_state

  END FUNCTION new_veg_memory

#endif
END MODULE mo_veg_memory
