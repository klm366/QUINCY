MODULE mo_jsb_memory
#ifndef __NO_JSBACH__

  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_util,              ONLY: one_of

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_jsbach_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory'

CONTAINS

  FUNCTION new_jsbach_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

    USE mo_jsb_varlist,       ONLY: t_var_list, jsb_get_var_list, jsb_new_var_list, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    ZAXIS_SURFACE, FILETYPE_NC2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid, new_vgrid
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model, get_vgrid, register_vgrid
    USE mo_jsb_base,          ONLY: get_model
    USE mo_atmland_memory,    ONLY: new_atm2land_memory, new_land2atm_memory
    USE mo_srf_memory,        ONLY: new_srf_memory
    USE mo_soil_memory,       ONLY: new_soil_memory
    USE mo_veg_memory,        ONLY: new_veg_memory
#ifndef __NO_JSBACH_HD__
    USE mo_hd_memory,         ONLY: new_hd_memory
#endif

    CHARACTER(len=*),             INTENT(in) :: varlist_name
    CHARACTER(len=*),             INTENT(in) :: prefix        ! Prefix for variable names
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

    CHARACTER(len=10) :: parent_tile_name

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_jsbach_memory'

    model => get_model(tile%owner_model_id)

    CALL message(TRIM(routine), 'Starting construction of new memory for '//TRIM(model%name)//' '//TRIM(tile%name))

    ALLOCATE(new)

    table = tables(1)

    CALL jsb_get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL jsb_new_var_list(this_list, TRIM(varlist_name),              &
                            loutput=.TRUE., output_type=FILETYPE_NC2,   &
                            lrestart=.TRUE., restart_type=FILETYPE_NC2, &
                            patch_id=model%grid%patch%id, table=table)
    END IF

    hgrid        => model%grid

    surface => get_vgrid(model, 'surface')
    IF (.NOT. ASSOCIATED(surface)) THEN
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(model, surface)
    END IF

    IF (ASSOCIATED(tile%parent_tile)) THEN
       parent_tile_name = tile%parent_tile%name
    ELSE
       parent_tile_name = 'grid box  '
    END IF

    CALL jsb_add_var(this_list, prefix//'fract'//suffix, new%fract,                                           &
      hgrid, surface,                                                                                         &
      t_cf('land_area_fraction', '',                                                                          &
               'fraction of '//TRIM(parent_tile_name)//' covered by '//TRIM(tile%name), nc_bits),             &
      t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      initval_r=0.0_wp)

    CALL jsb_add_var(this_list, prefix//'fract_max'//suffix, new%fract_max,                                   &
      hgrid, surface,                                                                                         &
      t_cf('land_area_fraction', '',                                                                          &
               'fraction of '//TRIM(parent_tile_name)//' covered by '//TRIM(tile%name), nc_bits),             &
      t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
      initval_r=0.0_wp)

    IF (ASSOCIATED(tile%parent_tile)) THEN
      ! For child tiles, atm2land just points to the memory of the root tile
      new%atm2land => tile%parent_tile%mem%atm2land
    ELSE
!!$      new%atm2land => new_atm2land_memory(varlist_name, 'a2l_'//TRIM(prefix), suffix, tile)
      new%atm2land => new_atm2land_memory(varlist_name, TRIM(prefix), suffix, tile)
!!$      new%land2atm => new_land2atm_memory(varlist_name, 'l2a_', suffix, tile)  ! land2atm is only necessary for the root tile
    END IF

    IF (model%config%srf%active .AND. one_of('srf', tile%components) > 0) &
!!$      new%srf => new_srf_memory(varlist_name, 'srf_'//TRIM(prefix), suffix, tile)
      new%srf => new_srf_memory(varlist_name, TRIM(prefix), suffix, tile)

    IF (model%config%soil%active .AND. one_of('soil', tile%components) > 0) &
!!$      new%soil => new_soil_memory(varlist_name, 'soil_'//TRIM(prefix), suffix, tile)
      new%soil => new_soil_memory(varlist_name, TRIM(prefix), suffix, tile)

    IF (model%config%veg%active .AND. one_of('veg', tile%components) > 0) &
!!$      new%veg => new_veg_memory(varlist_name, 'veg_'//TRIM(prefix), suffix, tile)
      new%veg => new_veg_memory(varlist_name, TRIM(prefix), suffix, tile)

#ifndef __NO_JSBACH_HD__
    IF (model%config%hd%active) &
      new%hd => new_hd_memory(varlist_name, TRIM(prefix), suffix, tile)
#endif

    mem => new

  END FUNCTION new_jsbach_memory

#endif
END MODULE mo_jsb_memory
