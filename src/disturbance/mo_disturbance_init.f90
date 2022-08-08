!>
!! @brief Initialization of the the disturbance memory.
!!
!! @author
!!  Rackhun Son, MPI-BGC Jena
!!
!! @par Revision History
!! Rackhun Son (2022-08-03): ligtnining climatology & population density for fire
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
      & lightning_climatology(:,:,:), &
      & population_density(:,:,:)

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
      & disturb_init_vars%population_density           &
      & )

  END SUBROUTINE disturb_finalize_init_vars

  SUBROUTINE disturb_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    dsl4jsb_Def_config(DISTURB_)

    REAL(wp), POINTER ::  &
      & ptr_3D(:,:,:)
    TYPE(t_input_file) :: input_file

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ALLOCATE( &
      & disturb_init_vars%lightning_climatology(nproma, nblks, 1:366),              &
      & disturb_init_vars%population_density   (nproma, nblks, 1:366)               &
      & )

    disturb_init_vars%lightning_climatology(:,:,:) = missval
    disturb_init_vars%population_density   (:,:,:) = missval

    dsl4jsb_Get_config(DISTURB_)

    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (lightning climatology) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%bc_lightning_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%bc_lightning_filename), model%grid_id)

    !  lightning climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading lightning_climatology ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='lightning_climatology')
    disturb_init_vars%lightning_climatology(:,:,1:366) = ptr_3D(:,:,1:366)
    DEALLOCATE(ptr_3D)

    CALL input_file%Close()

    IF (debug_on()) CALL message(TRIM(routine), 'Reading/setting disturbance (population density) init from ' &
      &                          //TRIM(dsl4jsb_Config(DISTURB_)%bc_population_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%bc_population_filename), model%grid_id)

    ! population density
    IF (debug_on()) CALL message(TRIM(routine), 'reading population density ...')
    ptr_3D => input_file%Read_2d_time(    &
      & variable_name='population_density')
    disturb_init_vars%population_density(:,:,1:366) = ptr_3D(:,:,1:366)
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
    !$acc enter data attach(dsl4jsb_memory(DISTURB_)%lightning_climatology)
    !$acc update device(dsl4jsb_memory(DISTURB_)%lightning_climatology)

    ! population density
    IF (debug_on()) CALL message(TRIM(routine), 'setting population density ...')
    dsl4jsb_memory(DISTURB_)%population_density(:,:,:) = disturb_init_vars%population_density(:,:,:)
    !$acc update device(dsl4jsb_memory(DISTURB_)%population_density)

  END SUBROUTINE disturb_init_bc

#endif
END MODULE mo_disturb_init


