!>
!! @brief Contains methods for the soil memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_memory

  USE mo_soil_memory_class, ONLY: t_soil_memory, t_soil_phys_state, &
                                  t_soil_energy_state, t_soil_energy_flux,          &
                                  t_soil_water_state, t_soil_water_flux
  USE mo_jsb_grid_class,    ONLY: get_nproma, get_nblks
  USE mo_cdi_constants,     ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16,  &
                                  ZA_SURFACE, ZA_DEPTH_BELOW_LAND,                  &
                                  FILETYPE_NC2,                                     &
                                  GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL,&
                                  TSTEP_CONSTANT
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_soil_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_memory'

CONTAINS

  FUNCTION new_soil_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

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
    TYPE(t_soil_memory), POINTER           :: mem

    TYPE(t_soil_memory), POINTER :: new
    TYPE(t_var_list),   POINTER :: this_list
    TYPE(t_jsb_model), POINTER  :: model

    INTEGER :: nproma, nblks
    INTEGER :: nsoil_energy, nsoil_water
    INTEGER :: shape2d(2), shape3d_e(3), shape3d_w(3)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_soil_memory'

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
    nsoil_energy = model%config%soil%nsoil_energy
    nsoil_water  = model%config%soil%nsoil_water

    shape2d   = (/ nproma, nblks /)
    shape3d_e = (/ nproma, nsoil_energy,nblks /)
    shape3d_w = (/ nproma, nsoil_water, nblks /)

    IF (ASSOCIATED(tile%parent_tile)) THEN
       new%phys_state => tile%parent_tile%mem%soil%phys_state
    ELSE
       new%phys_state => new_phys_state()
    END IF

    new%energy_state => new_energy_state()
    new%energy_flux  => new_energy_flux()
    new%water_state  => new_water_state()
    new%water_flux   => new_water_flux()

    mem => new

  CONTAINS

    FUNCTION new_phys_state() RESULT(pstate)

      TYPE(t_soil_phys_state), POINTER :: pstate

      ALLOCATE(pstate)

      CALL add_var( this_list, prefix//'soil_depth'//suffix, pstate%soil_depth,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('soil_depth', 'm', 'Depth of soil until bedrock', DATATYPE_FLT32),                   &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'soil_depth_l'//suffix, pstate%soil_depth_l,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                                  &
        & t_cf_var('soil_depth_layer', 'm', '', DATATYPE_FLT32),                                        &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape3d_w, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'soil_por'//suffix, pstate%vpor,                                 &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('vol_soil_porosity', '', '', DATATYPE_FLT32),                                        &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'pore_size'//suffix, pstate%psdi,                                &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('pore_size_dist_index', '', '', DATATYPE_FLT32),                                     &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'bclapp'//suffix, pstate%bclapp,                                 &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('exp_clapp_hornberger', '', '', DATATYPE_FLT32),                                     &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'matrix_pot'//suffix, pstate%mpot,                               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('matrix_potential', 'm', '', DATATYPE_FLT32),                                        &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'sat_hyd_cond'//suffix, pstate%shc,                              &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('sat_hydraulic_cond', 'm s-1', '', DATATYPE_FLT32),                                  &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'vol_heat_cap'//suffix, pstate%vhcap,                            &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('vol_heat_capacity_soil', 'J m-3 K-1', '', DATATYPE_FLT32),                          &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'heat_diff'//suffix, pstate%hdiff,                               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('heat_diffusivity_soil', 'm2 s-1', '', DATATYPE_FLT32),                              &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'alb_vis_soil'//suffix, pstate%alb_vis_soil,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('soil_albedo_visible ', '', '', DATATYPE_FLT32),                                     &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'alb_nir_soil'//suffix, pstate%alb_nir_soil,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('soil_albedo_nir', '', '', DATATYPE_FLT32),                                          &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_phys_state

    FUNCTION new_energy_state() RESULT(estate)

      TYPE(t_soil_energy_state), POINTER :: estate

      ALLOCATE(estate)

      CALL add_var( this_list, prefix//'t_soil'//suffix, estate%t_soil,                                 &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                                  &
        & t_cf_var('soil_temperature', 'K', '', DATATYPE_FLT32),                                        &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape3d_e, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'t_soil_aceof'//suffix, estate%t_soil_acoef,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                                  &
        & t_cf_var('soil_temperature_acoef', '', '', DATATYPE_FLT32),                                   &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape3d_e, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'t_soil_bceof'//suffix, estate%t_soil_bcoef,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                                  &
        & t_cf_var('soil_temperature_bcoef', '', '', DATATYPE_FLT32),                                   &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape3d_e, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'heat_cap_grnd'//suffix, estate%hcap_grnd,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('heat_capacity_ground', 'J m-2 K-1', '', DATATYPE_FLT32),                            &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_energy_state

    FUNCTION new_energy_flux() RESULT(eflux)

      TYPE(t_soil_energy_flux), POINTER :: eflux

      ALLOCATE(eflux)

      CALL add_var( this_list, prefix//'grnd_hflx'//suffix, eflux%grnd_hflx,                            &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('grnd_hflx', '', 'Ground heat flux', DATATYPE_FLT32),                                &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_energy_flux

    FUNCTION new_water_state() RESULT(wstate)

      TYPE(t_soil_water_state), POINTER :: wstate

      ALLOCATE(wstate)

      CALL add_var( this_list, prefix//'wfract_soil'//suffix, wstate%wfract_soil,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('wet_fraction_soil', '', '', DATATYPE_FLT32),                                        &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'wsr_soil'//suffix, wstate%wsr_soil,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('skin_reservoir_soil', 'm', '', DATATYPE_FLT32),                                     &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'sfract_soil'//suffix, wstate%sfract_soil,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('snow_fraction_soil', '', '', DATATYPE_FLT32),                                       &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'wsn_soil'//suffix, wstate%wsn_soil,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('snow_water_equivalent_soil', 'm', '', DATATYPE_FLT32),                              &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'ws_l'//suffix, wstate%ws_l,                                     &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                                  &
        & t_cf_var('soil_moisture', 'm', '', DATATYPE_FLT32),                                           &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape3d_w, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'ws_fc_l'//suffix, wstate%ws_fc_l,                               &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                                  &
        & t_cf_var('field_capacity', 'm', '', DATATYPE_FLT32),                                          &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape3d_w, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'ws_pwp_l'//suffix, wstate%ws_pwp_l,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                                  &
        & t_cf_var('permanent_wilting_point', 'm', '', DATATYPE_FLT32),                                 &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape3d_w, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'ws_sat_l'//suffix, wstate%ws_sat_l,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                                                  &
        & t_cf_var('water holding capacity', 'm', '', DATATYPE_FLT32),                                  &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape3d_w, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'ws'//suffix, wstate%ws,                                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('soil_moisture', 'm', '', DATATYPE_FLT32),                                           &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_water_state

    FUNCTION new_water_flux() RESULT(wflux)

      TYPE(t_soil_water_flux), POINTER :: wflux

      ALLOCATE(wflux)

      CALL add_var( this_list, prefix//'runoff'//suffix, wflux%runoff,                                  &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_runoff', '', '', DATATYPE_FLT32),                                           &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'drainage'//suffix, wflux%drainage,                              &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('drainage', '', '', DATATYPE_FLT32),                                                 &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_water_flux

  END FUNCTION new_soil_memory

END MODULE mo_soil_memory
