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
#ifndef __NO_JSBACH__

  USE mo_soil_memory_class, ONLY: t_soil_memory, t_soil_phys_state, &
                                  t_soil_energy_state, t_soil_energy_flux,          &
                                  t_soil_water_state, t_soil_water_flux
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_soil_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_memory'

CONTAINS

  FUNCTION new_soil_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

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
    TYPE(t_soil_memory), POINTER           :: mem

    TYPE(t_soil_memory), POINTER :: new
    TYPE(t_var_list),   POINTER :: this_list
    TYPE(t_jsb_model), POINTER  :: model

    INTEGER :: nsoil_energy, nsoil_water
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_e, soil_w      ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_soil_memory'

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

    nsoil_energy = model%config%soil%nsoil_energy
    nsoil_water  = model%config%soil%nsoil_water
    surface => get_vgrid(model, 'surface')
    IF (.NOT. ASSOCIATED(surface)) THEN
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(model, surface)
    END IF
    soil_e => get_vgrid(model, 'soil_depth_energy')
    IF (.NOT. ASSOCIATED(soil_e)) THEN
      soil_e  => new_vgrid('soil_depth_energy', ZAXIS_DEPTH_BELOW_LAND, nsoil_energy, &
                           levels=model%config%soil%soil_depth_energy(1:nsoil_energy), units='m')
      CALL register_vgrid(model, soil_e)
    END IF
    soil_w => get_vgrid(model, 'soil_depth_water')
    IF (.NOT. ASSOCIATED(soil_w)) THEN
      soil_w  => new_vgrid('soil_depth_water', ZAXIS_DEPTH_BELOW_LAND, nsoil_water, &
                           levels=model%config%soil%soil_depth_water(1:nsoil_water), units='m')
      CALL register_vgrid(model, soil_w)
    END IF

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

      CALL jsb_add_var( this_list, prefix//'soil_depth'//suffix, pstate%soil_depth,                         &
        & hgrid, surface,                                                                    &
        & t_cf('soil_depth', 'm', 'Depth of soil until bedrock', nc_bits),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'soil_depth_l'//suffix, pstate%soil_depth_l,                     &
        & hgrid, soil_w,                                                                     &
        & t_cf('soil_depth_layer', 'm', '', nc_bits),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'soil_por'//suffix, pstate%vpor,                                 &
        & hgrid, surface,                                                                    &
        & t_cf('vol_soil_porosity', '', '', nc_bits),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'pore_size'//suffix, pstate%psdi,                                &
        & hgrid, surface,                                                                    &
        & t_cf('pore_size_dist_index', '', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'bclapp'//suffix, pstate%bclapp,                                 &
        & hgrid, surface,                                                                    &
        & t_cf('exp_clapp_hornberger', '', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'matrix_pot'//suffix, pstate%mpot,                               &
        & hgrid, surface,                                                                    &
        & t_cf('matrix_potential', 'm', '', nc_bits),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'sat_hyd_cond'//suffix, pstate%shc,                              &
        & hgrid, surface,                                                                    &
        & t_cf('sat_hydraulic_cond', 'm s-1', '', nc_bits),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'vol_heat_cap'//suffix, pstate%vhcap,                            &
        & hgrid, surface,                                                                    &
        & t_cf('vol_heat_capacity_soil', 'J m-3 K-1', '', nc_bits),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'heat_diff'//suffix, pstate%hdiff,                               &
        & hgrid, surface,                                                                    &
        & t_cf('heat_diffusivity_soil', 'm2 s-1', '', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'alb_vis_soil'//suffix, pstate%alb_vis_soil,                     &
        & hgrid, surface,                                                                    &
        & t_cf('soil_albedo_visible ', '', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'alb_nir_soil'//suffix, pstate%alb_nir_soil,                     &
        & hgrid, surface,                                                                    &
        & t_cf('soil_albedo_nir', '', '', nc_bits),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

    END FUNCTION new_phys_state

    FUNCTION new_energy_state() RESULT(estate)

      TYPE(t_soil_energy_state), POINTER :: estate

      ALLOCATE(estate)

      CALL jsb_add_var( this_list, prefix//'t_soil'//suffix, estate%t_soil,                                 &
        & hgrid, soil_e,                                                                     &
        & t_cf('soil_temperature', 'K', '', nc_bits),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'t_soil_aceof'//suffix, estate%t_soil_acoef,                     &
        & hgrid, soil_e,                                                                     &
        & t_cf('soil_temperature_acoef', '', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'t_soil_bceof'//suffix, estate%t_soil_bcoef,                     &
        & hgrid, soil_e,                                                                     &
        & t_cf('soil_temperature_bcoef', '', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'heat_cap_grnd'//suffix, estate%hcap_grnd,                       &
        & hgrid, surface,                                                                    &
        & t_cf('heat_capacity_ground', 'J m-2 K-1', 'Ground heat capacity', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'heat_cap_grnd_old'//suffix, estate%hcap_grnd_old,&
        & hgrid, surface,                                                                    &
        & t_cf('heat_capacity_ground_old', 'J m-2 K-1', 'Ground heat capacity (old)', nc_bits), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

    END FUNCTION new_energy_state

    FUNCTION new_energy_flux() RESULT(eflux)

      TYPE(t_soil_energy_flux), POINTER :: eflux

      ALLOCATE(eflux)

      CALL jsb_add_var( this_list, prefix//'grnd_hflx'//suffix, eflux%grnd_hflx,                            &
        & hgrid, surface,                                                                    &
        & t_cf('grnd_hflx', '', 'Ground heat flux', nc_bits),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'grnd_hflx_old'//suffix, eflux%grnd_hflx_old,                    &
        & hgrid, surface,                                                                    &
        & t_cf('grnd_hflx_old', '', 'Ground heat flux (old)', nc_bits),                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

    END FUNCTION new_energy_flux

    FUNCTION new_water_state() RESULT(wstate)

      TYPE(t_soil_water_state), POINTER :: wstate

      ALLOCATE(wstate)

      CALL jsb_add_var( this_list, prefix//'wfract_soil'//suffix, wstate%wfract_soil,                       &
        & hgrid, surface,                                                                    &
        & t_cf('wet_fraction_soil', '', '', nc_bits),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsr_soil'//suffix, wstate%wsr_soil,                             &
        & hgrid, surface,                                                                    &
        & t_cf('skin_reservoir_soil', 'm', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'sfract_soil'//suffix, wstate%sfract_soil,                       &
        & hgrid, surface,                                                                    &
        & t_cf('snow_fraction_soil', '', '', nc_bits),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsn_soil'//suffix, wstate%wsn_soil,                             &
        & hgrid, surface,                                                                    &
        & t_cf('snow_water_equivalent_soil', 'm', '', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'overflow'//suffix, wstate%overflow,                             &
        & hgrid, surface,                                                                    &
        & t_cf('overflown_soil_moisture', 'm', '', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_l'//suffix, wstate%ws_l,                                     &
        & hgrid, soil_w,                                                                     &
        & t_cf('soil_moisture', 'm', '', nc_bits),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_fc_l'//suffix, wstate%ws_fc_l,                               &
        & hgrid, soil_w,                                                                     &
        & t_cf('field_capacity', 'm', '', nc_bits),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_pwp_l'//suffix, wstate%ws_pwp_l,                             &
        & hgrid, soil_w,                                                                     &
        & t_cf('permanent_wilting_point', 'm', '', nc_bits),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_sat_l'//suffix, wstate%ws_sat_l,                             &
        & hgrid, soil_w,                                                                     &
        & t_cf('water holding capacity', 'm', '', nc_bits),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws'//suffix, wstate%ws,                                         &
        & hgrid, surface,                                                                    &
        & t_cf('soil_moisture', 'm', '', nc_bits),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

    END FUNCTION new_water_state

    FUNCTION new_water_flux() RESULT(wflux)

      TYPE(t_soil_water_flux), POINTER :: wflux

      ALLOCATE(wflux)

      CALL jsb_add_var( this_list, prefix//'runoff'//suffix, wflux%runoff,                   &
        & hgrid, surface,                                                                    &
        & t_cf('surface_runoff', '', '', nc_bits),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'drainage'//suffix, wflux%drainage,               &
        & hgrid, surface,                                                                    &
        & t_cf('drainage', '', '', nc_bits),                                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

    END FUNCTION new_water_flux

  END FUNCTION new_soil_memory

#endif
END MODULE mo_soil_memory
