!>
!! @brief Contains methods for the srf memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-27)
!!
MODULE mo_srf_memory

  USE mo_srf_memory_class, ONLY: t_srf_memory, t_srf_phys_state, &
                                  t_srf_energy_state, t_srf_energy_flux,          &
                                  t_srf_water_state,  t_srf_water_flux
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

  PUBLIC :: new_srf_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_memory'

CONTAINS

  FUNCTION new_srf_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

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
    TYPE(t_srf_memory), POINTER           :: mem

    TYPE(t_srf_memory), POINTER :: new
    TYPE(t_var_list),   POINTER :: this_list
    TYPE(t_jsb_model), POINTER  :: model

    INTEGER :: nproma, nblks
    INTEGER :: shape2d(2)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_srf_memory'

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

    shape2d   = (/ nproma, nblks /)

    IF (ASSOCIATED(tile%parent_tile)) THEN
       new%phys_state => tile%parent_tile%mem%srf%phys_state
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

      TYPE(t_srf_phys_state), POINTER :: pstate

      ALLOCATE(pstate)

      CALL add_var( this_list, prefix//'oro_stddev'//suffix, pstate%oro_stddev,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('orography_standard_deviation', '', '', DATATYPE_FLT32),                             &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'albedo_srf'//suffix, pstate%albedo_srf,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_albedo_background', '', '', DATATYPE_FLT32),                                &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL add_var( this_list, prefix//'alb_vis_srf'//suffix, pstate%alb_vis_srf,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_albedo_visible', '', '', DATATYPE_FLT32),                                   &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'alb_nir_srf'//suffix, pstate%alb_nir_srf,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_albedo_nir', '', '', DATATYPE_FLT32),                                       &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_phys_state

    FUNCTION new_energy_state() RESULT(estate)

      TYPE(t_srf_energy_state), POINTER :: estate

      ALLOCATE(estate)

      CALL add_var( this_list, prefix//'t_srf'//suffix, estate%t_srf,                                   &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_temperature', 'K', '', DATATYPE_FLT32),                                     &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

!!$      CALL add_var( this_list, prefix//'t_srf_old'//suffix, estate%t_srf_old,                           &
!!$        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
!!$        & t_cf_var('surface_temperature_old', 'K', '', DATATYPE_FLT32),                                 &
!!$        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
!!$        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'t_srf_unfilt'//suffix, estate%t_srf_unfilt,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_temperature_unfiltered', 'K', '', DATATYPE_FLT32),                          &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'t_rad_srf'//suffix, estate%t_rad_srf,                           &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_temperature_rad', 'K', '', DATATYPE_FLT32),                                 &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

!!$      CALL add_var( this_list, prefix//'t_eff_srf'//suffix, estate%t_eff_srf,                           &
!!$        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
!!$        & t_cf_var('surface_temperature_eff', 'K', '', DATATYPE_FLT32),                                 &
!!$        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
!!$        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'s_srf'//suffix, estate%s_srf,                                   &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('dry_static_energy', '', '', DATATYPE_FLT32),                                        &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

!!$      CALL add_var( this_list, prefix//'s_srf_old'//suffix, estate%s_srf_old,                           &
!!$        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
!!$        & t_cf_var('dry_static_energy_old', '', '', DATATYPE_FLT32),                                    &
!!$        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
!!$        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'qsat_srf'//suffix, estate%qsat_srf,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('sat_specific_humidity', '', '', DATATYPE_FLT32),                                    &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

!!$      CALL add_var( this_list, prefix//'qsat_srf_old'//suffix, estate%qsat_srf_old,                     &
!!$        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
!!$        & t_cf_var('sat_specific_humidity_old', '', '', DATATYPE_FLT32),                                &
!!$        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
!!$        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'fact_q_air'//suffix, estate%fact_q_air,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('factor_air_moisture', '', '', DATATYPE_FLT32),                                      &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'fact_qsat_srf'//suffix, estate%fact_qsat_srf,                   &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('factor_qsat_srf', '', '', DATATYPE_FLT32),                                          &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_energy_state

    FUNCTION new_energy_flux() RESULT(eflux)

      TYPE(t_srf_energy_flux), POINTER :: eflux

      ALLOCATE(eflux)

      CALL add_var( this_list, prefix//'swrad_srf_net'//suffix, eflux%swrad_srf_net,                    &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('swrad_srf_net', '', 'Net surface shortwave radiation', DATATYPE_FLT32),             &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'lwrad_srf_net'//suffix, eflux%lwrad_srf_net,                    &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('lwrad_srf_net', '', 'Net surface longwave radiation', DATATYPE_FLT32),              &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'rad_srf_net'//suffix, eflux%rad_srf_net,                        &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('rad_srf_net', '', 'Net surface radiation', DATATYPE_FLT32),                         &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'latent_hflx'//suffix, eflux%latent_hflx,                        &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('latent_hflx', '', 'Latent heat flux', DATATYPE_FLT32),                              &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_energy_flux

    FUNCTION new_water_state() RESULT(wstate)

      TYPE(t_srf_water_state), POINTER :: wstate

      ALLOCATE(wstate)

      CALL add_var( this_list, prefix//'wfract_srf'//suffix, wstate%wfract_srf,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_wet_fraction', '', '', DATATYPE_FLT32),                                     &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'wsr_srf'//suffix, wstate%wsr_srf,                               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_skin_reservoir', 'm', '', DATATYPE_FLT32),                                  &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'sfract_srf'//suffix, wstate%sfract_srf,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_snow_fraction', '', '', DATATYPE_FLT32),                                    &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

      CALL add_var( this_list, prefix//'wsn_srf'//suffix, wstate%wsn_srf,                               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                           &
        & t_cf_var('surface_snow_water_equivalent', 'm', '', DATATYPE_FLT32),                           &
        & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),                       &
        & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_water_state

    FUNCTION new_water_flux() RESULT(wflux)

      TYPE(t_srf_water_flux), POINTER :: wflux

      ALLOCATE(wflux)

      CALL add_var( this_list, prefix//'evap_srf'//suffix, wflux%evap_srf,                              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                                 &
                  & t_cf_var('surface_evaporation', '', 'Evaporation from surface', DATATYPE_FLT32),    &
                  & t_grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_REFERENCE, GRID_CELL),             &
                  & ldims=shape2d, initval_r=0.0_wp )

    END FUNCTION new_water_flux

  END FUNCTION new_srf_memory

END MODULE mo_srf_memory
