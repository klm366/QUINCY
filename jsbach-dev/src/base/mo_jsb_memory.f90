MODULE mo_jsb_memory

  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_cdi_constants,     ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, &
       ZA_SURFACE, &
       FILETYPE_NC2, &
       GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, print_value
  USE mo_util_string,       ONLY: one_of

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_jsbach_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory'

CONTAINS

  FUNCTION new_jsbach_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

    USE mo_jsb_varlist_iface, ONLY: t_var_list, get_var_list, new_var_list, default_var_list_settings, add_var
    USE mo_jsb_io_iface,      ONLY: t_cf_var, t_grib2_var
    USE mo_jsb_base,          ONLY: get_model
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, get_nproma, get_nblks
    USE mo_atmland_memory,    ONLY: new_atm2land_memory, new_land2atm_memory
    USE mo_srf_memory,        ONLY: new_srf_memory
    USE mo_soil_memory,       ONLY: new_soil_memory
    USE mo_veg_memory,        ONLY: new_veg_memory

    CHARACTER(len=*),             INTENT(in) :: varlist_name
    CHARACTER(len=*),             INTENT(in) :: prefix
    CHARACTER(len=*),             INTENT(in) :: suffix
#ifdef HAVE_F2003
    TYPE(t_jsb_tile),   POINTER,  INTENT(in) :: tile
#else
    TYPE(t_jsb_tile),   POINTER              :: tile
#endif
    TYPE(t_jsb_memory), POINTER              :: mem

    TYPE(t_jsb_memory), POINTER :: new
    TYPE(t_var_list),   POINTER :: this_list
    TYPE(t_jsb_model),  POINTER :: model

    INTEGER :: shape2d(2)
    CHARACTER(len=10) :: parent_tile_name

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_jsbach_memory'

    model => get_model(tile%owner_model_id)

    CALL message(TRIM(routine), 'Starting construction of new memory for '//TRIM(model%name)//' '//TRIM(tile%name))

    ALLOCATE(new)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL new_var_list(this_list, TRIM(varlist_name), patch_id=model%grid%patch%id)
      CALL default_var_list_settings(this_list, lrestart=.TRUE., restart_type=FILETYPE_NC2)
    END IF

    shape2d = (/ get_nproma(model%grid), get_nblks(model%grid) /)
    IF (ASSOCIATED(tile%parent_tile)) THEN
       parent_tile_name = tile%parent_tile%name
    ELSE
       parent_tile_name = 'grid box  '
    END IF

    CALL add_var(this_list, prefix//'fract'//suffix, new%fract,                                               &
      GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                                     &
      t_cf_var('land_area_fraction', '',                                                                      &
               'fraction of '//TRIM(parent_tile_name)//' covered by '//TRIM(tile%name), DATATYPE_FLT32),      &
      t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                                 &
      ldims=shape2d, initval_r=0.0_wp)

    CALL add_var(this_list, prefix//'fract_max'//suffix, new%fract_max,                                       &
      GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                                     &
      t_cf_var('land_area_fraction', '',                                                                      &
               'fraction of '//TRIM(parent_tile_name)//' covered by '//TRIM(tile%name), DATATYPE_FLT32),      &
      t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                                 &
      ldims=shape2d, initval_r=0.0_wp)

    IF (ASSOCIATED(tile%parent_tile)) THEN
      ! For child tiles, atm2land just points to the memory of the root tile
      new%atm2land => tile%parent_tile%mem%atm2land
    ELSE
!!$      NEW%atm2land => new_atm2land_memory(varlist_name, 'a2l_'//TRIM(prefix), suffix, tile)
      NEW%atm2land => new_atm2land_memory(varlist_name, TRIM(prefix), suffix, tile)
!!$      new%land2atm => new_land2atm_memory(varlist_name, 'l2a_', suffix, tile)  ! land2atm is only necessary for the root tile
    END IF

    IF (model%config%srf%active .AND. one_of('srf', tile%components) > 0) &
!!$      NEW%srf => new_srf_memory(varlist_name, 'srf_'//TRIM(prefix), suffix, tile)
      NEW%srf => new_srf_memory(varlist_name, TRIM(prefix), suffix, tile)

    IF (model%config%soil%active .AND. one_of('soil', tile%components) > 0) &
!!$      NEW%soil => new_soil_memory(varlist_name, 'soil_'//TRIM(prefix), suffix, tile)
      NEW%soil => new_soil_memory(varlist_name, TRIM(prefix), suffix, tile)

    IF (model%config%veg%active .AND. one_of('veg', tile%components) > 0) &
!!$      NEW%veg => new_veg_memory(varlist_name, 'veg_'//TRIM(prefix), suffix, tile)
      NEW%veg => new_veg_memory(varlist_name, TRIM(prefix), suffix, tile)

    mem => new

  END FUNCTION new_jsbach_memory

END MODULE mo_jsb_memory
