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
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_veg_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory'

CONTAINS

  FUNCTION new_veg_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

    USE mo_jsb_varlist,       ONLY: t_var_list, jsb_get_var_list, jsb_new_var_list, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    FILETYPE_NC2, TSTEP_CONSTANT, ZAXIS_SURFACE, ZAXIS_DEPTH_BELOW_LAND, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid, new_vgrid
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model, get_vgrid, register_vgrid
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

    INTEGER                    :: nsoil_water
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_w              ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_veg_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(new)

    table = tables(1)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL jsb_get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL jsb_new_var_list(this_list, TRIM(varlist_name),              &
                            lrestart=.TRUE., restart_type=FILETYPE_NC2, &
                            patch_id=model%grid%patch%id, table=table)
    END IF

    hgrid        => model%grid

    surface => get_vgrid(model, 'surface')
    IF (.NOT. ASSOCIATED(surface)) THEN
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(model, surface)
    END IF
    soil_w => get_vgrid(model, 'soil_depth_water')
    IF (.NOT. ASSOCIATED(soil_w)) THEN
      nsoil_water  = model%config%soil%nsoil_water
      soil_w  => new_vgrid('soil_depth_water', ZAXIS_DEPTH_BELOW_LAND, nsoil_water, &
                           levels=model%config%soil%soil_depth_water(1:nsoil_water), units='m')
      CALL register_vgrid(model, soil_w)
    END IF

    new%phys_state  => new_phys_state()
    new%water_state => new_water_state()

    mem => new

  CONTAINS

    FUNCTION new_phys_state() RESULT(pstate)

      TYPE(t_veg_phys_state), POINTER :: pstate

      ALLOCATE(pstate)

      CALL jsb_add_var( this_list, prefix//'lai'//suffix, pstate%lai,                             &
        & hgrid, surface,                                                 &
        & t_cf('leaf_area_index', '', '', nc_bits),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      ! Note: extra dimensions (except levels) are not allowed in echam streams
!!$      CALL jsb_add_var( this_list, prefix//'lai_mon'//suffix, pstate%lai_mon,                     &
!!$        & hgrid, surface,                                                 &
!!$        & t_cf('monthly_leaf_area_index', '', '', nc_bits),                        &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
!!$        & extra_dim_n=14, initval_r=0.0_wp, loutput=.FALSE., lrestart=.FALSE. )
      ALLOCATE(pstate%lai_mon(hgrid%nproma, hgrid%nblks, 14))

      CALL jsb_add_var( this_list, prefix//'fract_fpc'//suffix, pstate%fract_fpc,                 &
        & hgrid, surface,                                                 &
        & t_cf('foliage_projected_cover', '', '', nc_bits),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'fract_fpc_max'//suffix, pstate%fract_fpc_max,         &
        & hgrid, surface,                                                 &
        & t_cf('maximum_fpc', '', '', nc_bits),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'fract_forest'//suffix, pstate%fract_forest,           &
        & hgrid, surface,                                                 &
        & t_cf('forest_fraction', '', '', nc_bits),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'canopy_cond'//suffix, pstate%canopy_cond,             &
        & hgrid, surface,                                                 &
        & t_cf('canopy_conductance', '', '', nc_bits),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'skin_can_depth'//suffix, pstate%skin_can_depth,       &
        & hgrid, surface,                                                 &
        & t_cf('canopy_skin_depth', 'm', '', nc_bits),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'root_depth'//suffix, pstate%root_depth,               &
        & hgrid, surface,                                                 &
        & t_cf('rooting_depth', '', '', nc_bits),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'root_depth_l'//suffix, pstate%root_depth_l,           &
        & hgrid, soil_w,                                                                      &
        & t_cf('root_depth_layer', 'm', '', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'alb_vis_can'//suffix, pstate%alb_vis_can,             &
        & hgrid, surface,                                                 &
        & t_cf('canopy_albedo_visible', '', '', nc_bits),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'alb_nir_can'//suffix, pstate%alb_nir_can,             &
        & hgrid, surface,                                                 &
        & t_cf('canopy_albedo_nir', '', '', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

    END FUNCTION new_phys_state

    FUNCTION new_water_state() RESULT(wstate)

      TYPE(t_veg_water_state), POINTER :: wstate

      ALLOCATE(wstate)

      CALL jsb_add_var( this_list, prefix//'wfract_can'//suffix, wstate%wfract_can,                         &
        & hgrid, surface,                                                           &
        & t_cf('wet_fraction_canopy', '', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsr_can'//suffix, wstate%wsr_can,                               &
        & hgrid, surface,                                                           &
        & t_cf('skin_reservoir_canopy', 'm', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'sfract_can'//suffix, wstate%sfract_can,                         &
        & hgrid, surface,                                                           &
        & t_cf('snow_fraction_canopy', '', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsn_can'//suffix, wstate%wsn_can,                               &
        & hgrid, surface,                                                           &
        & t_cf('snow_water_equivalent_canopy', 'm', '', nc_bits),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_root'//suffix, wstate%ws_root,                               &
        & hgrid, surface,                                                           &
        & t_cf('soil_moisture_root', 'm', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_fc_root'//suffix, wstate%ws_fc_root,                         &
        & hgrid, surface,                                                           &
        & t_cf('field_capacity_root', 'm', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_pwp_root'//suffix, wstate%ws_pwp_root,                       &
        & hgrid, surface,                                                           &
        & t_cf('wilting_root', 'm', '', nc_bits),                                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'water_stress'//suffix, wstate%water_stress,                     &
        & hgrid, surface,                                                           &
        & t_cf('water_stress_factor', '', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

    END FUNCTION new_water_state

  END FUNCTION new_veg_memory

#endif
END MODULE mo_veg_memory
