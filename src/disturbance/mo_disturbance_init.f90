!>
!! @brief Initialization of the the disturbance memory.
!!
!! @author
!!  Rackhun Son, MPI-BGC Jena
!!
!! @par Revision History
!! Rackhun Son (2022-08-03): ligtnining climatology & population density for fire
!! Rackhun Son (2022-09-11): local mean (climatology) variables for fire
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_disturb_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_JSBACH, MODEL_JSBACH_QUINCY, MODEL_QUINCY
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,              ONLY: missval

  dsl4jsb_Use_processes DISTURB_
  dsl4jsb_Use_config(DISTURB_)
  dsl4jsb_Use_memory(DISTURB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: disturb_init

  TYPE t_disturb_init_vars
    REAL(wp), POINTER ::  &
      & lightning_climatology(:,:,:),   &

      & slope(:,:),                     &
      & roughness(:,:),                 &
      & topographic_position_idx(:,:),  &
      & terrain_ruggedness_idx(:,:),    &
      & vector_ruggedness_measure(:,:), &
      & profile_curvature(:,:),         &
      & tangential_curvature(:,:),      &

      !& population_density(:,:,:), &

      & t_air_local_clim(:,:,:),        &
      & q_rel_local_clim(:,:,:),        &
      & q_local_clim(:,:,:),            &
      & ws10_local_clim(:,:,:),         &
      & swvl1_local_clim(:,:,:),        &
      & swvl2_local_clim(:,:,:),        &
      & swvl3_local_clim(:,:,:),        &
      & swvl4_local_clim(:,:,:),        &
      & lai_local_clim(:,:,:)

  END TYPE t_disturb_init_vars

  TYPE(t_disturb_init_vars) :: disturb_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_init'

CONTAINS

  !> Run disturbance init
  !! 
  SUBROUTINE disturb_init(tile)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local
    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_init'

    model => Get_model(tile%owner_model_id)

    ! currently DL based Fire model is only tested in JSBACH 
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_JSBACH)
      IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
        CALL disturb_read_init_vars(tile)
      END IF
      CALL disturb_init_bc(tile)
    END SELECT

    IF (tile%Is_last_process_tile(DISTURB_)) THEN
        CALL disturb_finalize_init_vars()
    END IF

  END SUBROUTINE disturb_init

  SUBROUTINE disturb_finalize_init_vars

    DEALLOCATE( &
      & disturb_init_vars%lightning_climatology,       &
      !& disturb_init_vars%population_density,          &

      & disturb_init_vars%slope,                       &
      & disturb_init_vars%roughness,                   &
      & disturb_init_vars%topographic_position_idx,    &
      & disturb_init_vars%terrain_ruggedness_idx,      &
      & disturb_init_vars%vector_ruggedness_measure,   &
      & disturb_init_vars%profile_curvature,           &
      & disturb_init_vars%tangential_curvature,        &

      & disturb_init_vars%t_air_local_clim,            &
      & disturb_init_vars%q_rel_local_clim,            &
      & disturb_init_vars%q_local_clim,                &
      & disturb_init_vars%ws10_local_clim,             &
      & disturb_init_vars%swvl1_local_clim,            &
      & disturb_init_vars%swvl2_local_clim,            &
      & disturb_init_vars%swvl3_local_clim,            &
      & disturb_init_vars%swvl4_local_clim,            &
      & disturb_init_vars%lai_local_clim               &

      & )

  END SUBROUTINE disturb_finalize_init_vars

  SUBROUTINE disturb_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    dsl4jsb_Def_config(DISTURB_)

    REAL(wp), POINTER ::  &
      & ptr_3D(:,:,:),    &     
      & ptr_2D(:,:)
    TYPE(t_input_file) :: input_file

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ALLOCATE( &
      & disturb_init_vars%lightning_climatology(nproma, nblks, 1:366),               &
      !& disturb_init_vars%population_density   (nproma, nblks, 1:366),               &

      & disturb_init_vars%slope(nproma, nblks),                                      &
      & disturb_init_vars%roughness(nproma, nblks),                                  &
      & disturb_init_vars%topographic_position_idx(nproma, nblks),                   &
      & disturb_init_vars%terrain_ruggedness_idx(nproma, nblks),                     &
      & disturb_init_vars%vector_ruggedness_measure(nproma, nblks),                  &
      & disturb_init_vars%profile_curvature(nproma, nblks),                          &
      & disturb_init_vars%tangential_curvature(nproma, nblks),                       &

      & disturb_init_vars%t_air_local_clim (nproma, nblks, 1:366),                   &
      & disturb_init_vars%q_rel_local_clim (nproma, nblks, 1:366),                   &
      & disturb_init_vars%q_local_clim     (nproma, nblks, 1:366),                   &
      & disturb_init_vars%ws10_local_clim  (nproma, nblks, 1:366),                   &
      & disturb_init_vars%swvl1_local_clim (nproma, nblks, 1:366),                   &
      & disturb_init_vars%swvl2_local_clim (nproma, nblks, 1:366),                   &
      & disturb_init_vars%swvl3_local_clim (nproma, nblks, 1:366),                   &
      & disturb_init_vars%swvl4_local_clim (nproma, nblks, 1:366),                   &
      & disturb_init_vars%lai_local_clim   (nproma, nblks, 1:366)                    &
      !& disturb_init_vars%t_air_local_clim (nproma, nblks, 1:12),                   &
      !& disturb_init_vars%q_rel_local_clim (nproma, nblks, 1:12),                   &
      !& disturb_init_vars%q_local_clim     (nproma, nblks, 1:12),                   &
      !& disturb_init_vars%ws10_local_clim  (nproma, nblks, 1:12),                   &
      !& disturb_init_vars%swvl1_local_clim (nproma, nblks, 1:12),                   &
      !& disturb_init_vars%swvl2_local_clim (nproma, nblks, 1:12),                   &
      !& disturb_init_vars%swvl3_local_clim (nproma, nblks, 1:12),                   &
      !& disturb_init_vars%swvl4_local_clim (nproma, nblks, 1:12),                   &
      !& disturb_init_vars%lai_local_clim   (nproma, nblks, 1:12)                    &
      & )

    disturb_init_vars%lightning_climatology(:,:,:) = missval
    !disturb_init_vars%population_density   (:,:,:) = missval

    disturb_init_vars%slope                    (:,:) = missval
    disturb_init_vars%roughness                (:,:) = missval
    disturb_init_vars%topographic_position_idx (:,:) = missval
    disturb_init_vars%terrain_ruggedness_idx   (:,:) = missval
    disturb_init_vars%vector_ruggedness_measure(:,:) = missval
    disturb_init_vars%profile_curvature        (:,:) = missval
    disturb_init_vars%tangential_curvature     (:,:) = missval

    disturb_init_vars%t_air_local_clim (:,:,:) = missval
    disturb_init_vars%q_rel_local_clim (:,:,:) = missval
    disturb_init_vars%q_local_clim     (:,:,:) = missval
    disturb_init_vars%ws10_local_clim  (:,:,:) = missval
    disturb_init_vars%swvl1_local_clim (:,:,:) = missval
    disturb_init_vars%swvl2_local_clim (:,:,:) = missval
    disturb_init_vars%swvl3_local_clim (:,:,:) = missval
    disturb_init_vars%swvl4_local_clim (:,:,:) = missval
    disturb_init_vars%lai_local_clim   (:,:,:) = missval

    dsl4jsb_Get_config(DISTURB_)

    !  lightning climatology
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (lightning climatology) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%bc_lightning_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%bc_lightning_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading lightning_climatology ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='lightning_climatology')
    disturb_init_vars%lightning_climatology(:,:,1:366) = ptr_3D(:,:,1:366)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    !  Topographical complexity
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (topographical complexity) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%bc_topography_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%bc_topography_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading topographical complexity variables ...')
    ptr_2D => input_file%Read_2d(    &
      & variable_name='slope')
    disturb_init_vars%slope(:,:) = ptr_2D(:,:)
    DEALLOCATE(ptr_2D)

    ptr_2D => input_file%Read_2d(    &
      & variable_name='roughness')
    disturb_init_vars%roughness(:,:) = ptr_2D(:,:)
    DEALLOCATE(ptr_2D)

    ptr_2D => input_file%Read_2d(    &
      & variable_name='tpi')
    disturb_init_vars%topographic_position_idx(:,:) = ptr_2D(:,:)
    DEALLOCATE(ptr_2D)

    ptr_2D => input_file%Read_2d(    &
      & variable_name='tri')
    disturb_init_vars%terrain_ruggedness_idx(:,:) = ptr_2D(:,:)
    DEALLOCATE(ptr_2D)

    ptr_2D => input_file%Read_2d(    &
      & variable_name='vrm')
    disturb_init_vars%vector_ruggedness_measure(:,:) = ptr_2D(:,:)
    DEALLOCATE(ptr_2D)

    ptr_2D => input_file%Read_2d(    &
      & variable_name='pcurv')
    disturb_init_vars%profile_curvature(:,:) = ptr_2D(:,:)
    DEALLOCATE(ptr_2D)

    ptr_2D => input_file%Read_2d(    &
      & variable_name='tcurv')
    disturb_init_vars%tangential_curvature(:,:) = ptr_2D(:,:)
    DEALLOCATE(ptr_2D)

    CALL input_file%Close()

    ! population density
    !IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (population density) init from ' &
    !  &                          //TRIM(dsl4jsb_Config(DISTURB_)%bc_population_filename))
    !
    !input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%bc_population_filename), model%grid_id)
    !
    !IF (debug_on()) CALL message(TRIM(routine), 'reading population density ...')
    !ptr_3D => input_file%Read_2d_time(    &
    !  & variable_name='population_density')
    !disturb_init_vars%population_density(:,:,1:366) = ptr_3D(:,:,1:366)
    !DEALLOCATE(ptr_3D)

    ! temperature (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (temperature local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_t_air_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_t_air_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean of temperature ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='t2m')
    disturb_init_vars%t_air_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%t_air_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    ! relative humidity (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (relative humidity local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_rH_cf_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_rH_cf_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean relative humidity ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='rH_cf')
    disturb_init_vars%q_rel_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%q_rel_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    ! specific humidity (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (specific humidity local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_q_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_q_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean specific humidity ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='q')
    disturb_init_vars%q_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%q_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    ! wind speed (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (wind speed local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_ws10_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_ws10_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean wind speed ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='ws10')
    disturb_init_vars%ws10_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%ws10_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    ! soil water content level1 (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (soil water level1 local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_swvl1_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_swvl1_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean soil water content level1 ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='swvl1')
    disturb_init_vars%swvl1_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%swvl1_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    ! soil water content level2 (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (soil water level2 local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_swvl2_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_swvl2_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean soil water content level2 ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='swvl2')
    disturb_init_vars%swvl2_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%swvl2_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    ! soil water content level3 (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (soil water level3 local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_swvl3_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_swvl3_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean soil water content level3 ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='swvl3')
    disturb_init_vars%swvl3_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%swvl3_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    ! soil water content level4 (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (soil water level4 local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_swvl4_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_swvl4_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean soil water content level4 ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='swvl4')
    disturb_init_vars%swvl4_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%swvl4_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    ! LAI (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (lai local mean) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_lai_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%dl_local_clim_lai_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading local mean lai ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='Lai')
    disturb_init_vars%lai_local_clim(:,:,1:366) = ptr_3D(:,:,1:366)
    !disturb_init_vars%lai_local_clim(:,:,1:12) = ptr_3D(:,:,1:12)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

  END SUBROUTINE disturb_read_init_vars

  SUBROUTINE disturb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(DISTURB_)
    dsl4jsb_Def_memory(DISTURB_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(DISTURB_)
    dsl4jsb_Get_memory(DISTURB_)

    IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of disturb memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(DISTURB_)%bc_lightning_filename))

    ! Lightning climatology
    IF (debug_on()) CALL message(TRIM(routine), 'setting lightning climatology ...')
    dsl4jsb_memory(DISTURB_)%lightning_climatology(:,:,:) = disturb_init_vars%lightning_climatology(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%lightning_climatology)

    ! Topographical complexity
    IF (debug_on()) CALL message(TRIM(routine), 'setting topographical complexities ...')

    dsl4jsb_memory(DISTURB_)%slope(:,:)                     = disturb_init_vars%slope(:,:)
    dsl4jsb_memory(DISTURB_)%roughness(:,:)                 = disturb_init_vars%roughness(:,:)
    dsl4jsb_memory(DISTURB_)%topographic_position_idx(:,:)  = disturb_init_vars%topographic_position_idx(:,:)
    dsl4jsb_memory(DISTURB_)%terrain_ruggedness_idx(:,:)    = disturb_init_vars%terrain_ruggedness_idx(:,:)
    dsl4jsb_memory(DISTURB_)%vector_ruggedness_measure(:,:) = disturb_init_vars%vector_ruggedness_measure(:,:)
    dsl4jsb_memory(DISTURB_)%profile_curvature(:,:)         = disturb_init_vars%profile_curvature(:,:)
    dsl4jsb_memory(DISTURB_)%tangential_curvature(:,:)      = disturb_init_vars%tangential_curvature(:,:)
    
    ! population density
    !IF (debug_on()) CALL message(TRIM(routine), 'setting population density ...')
    !dsl4jsb_memory(DISTURB_)%population_density(:,:,:) = disturb_init_vars%population_density(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%population_density)

    ! temperature (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting temperature local mean ...')
    dsl4jsb_memory(DISTURB_)%t_air_local_clim(:,:,:) = disturb_init_vars%t_air_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%t_air_local_clim)

    ! relative humidity (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting relative humidity local mean ...')
    dsl4jsb_memory(DISTURB_)%q_rel_local_clim(:,:,:) = disturb_init_vars%q_rel_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%q_rel_local_clim)

    ! speficif humidity (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting specific humidity local mean ...')
    dsl4jsb_memory(DISTURB_)%q_local_clim(:,:,:) = disturb_init_vars%q_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%q_local_clim)

    ! wind speed (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting wind speed local mean ...')
    dsl4jsb_memory(DISTURB_)%ws10_local_clim(:,:,:) = disturb_init_vars%ws10_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%ws10_local_clim)

    ! soil water content level1 (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting soil water content level1 local mean ...')
    dsl4jsb_memory(DISTURB_)%swvl1_local_clim(:,:,:) = disturb_init_vars%swvl1_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%swvl1_local_clim)

    ! soil water content level2 (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting soil water content level2 local mean ...')
    dsl4jsb_memory(DISTURB_)%swvl2_local_clim(:,:,:) = disturb_init_vars%swvl2_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%swvl2_local_clim)

    ! soil water content level3 (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting soil water content level3 local mean ...')
    dsl4jsb_memory(DISTURB_)%swvl3_local_clim(:,:,:) = disturb_init_vars%swvl3_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%swvl3_local_clim)

    ! soil water content level4 (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting soil water content level4 local mean ...')
    dsl4jsb_memory(DISTURB_)%swvl4_local_clim(:,:,:) = disturb_init_vars%swvl4_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%swvl4_local_clim)

    ! lai (local mean)
    IF (debug_on()) CALL message(TRIM(routine), 'setting lai local mean ...')
    dsl4jsb_memory(DISTURB_)%lai_local_clim(:,:,:) = disturb_init_vars%lai_local_clim(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%lai_local_clim)


  END SUBROUTINE disturb_init_bc

#endif
END MODULE mo_disturb_init


