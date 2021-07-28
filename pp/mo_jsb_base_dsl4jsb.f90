!
! This file was automatically generated from:
!    /scratch/snx3000/rschnur/icon-rc-testing-plus-mc8.test/externals/jsbach/src/base/mo_jsb_base.f90
!    with md5 hash: 58076fd46afb4578fa9f2839bbfc3f3c
!
!>
!! @brief Contains main JSBACH structure and methods.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_base
#ifndef __NO_JSBACH__

  USE mo_jsb_class,   ONLY: jsbach
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message
  USE mo_util,        ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach
  PUBLIC :: jsbach_setup_models, jsbach_setup_tiles

CONTAINS

  SUBROUTINE jsbach_setup_models(master_namelist_filename)

    USE mo_jsb_control,         ONLY: jsb_models_nml, get_no_of_models, init_jsb_master_control, &
      &                               jsbach_runs_standalone, jsbach_is_restarted
    USE mo_jsb_parallel,        ONLY: init_parallel_iface
    USE mo_jsb_model_class,     ONLY: t_jsb_model, new_model
    USE mo_jsb_config_class,    ONLY: new_model_config
    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process
    USE mo_jsb_time,            ONLY: init_time
    USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,            ONLY: Register_vgrid
    USE mo_jsb_io,              ONLY: ZAXIS_SURFACE

    USE mo_jsb_process_class, ONLY: sse_, hydro_, rad_, carbon_, disturb_, fuel_
    USE mo_sse_config_class, ONLY: t_sse_config
    USE mo_hydro_config_class, ONLY: t_hydro_config
    USE mo_rad_config_class, ONLY: t_rad_config
    USE mo_carbon_config_class, ONLY: t_carbon_config
    USE mo_disturb_config_class, ONLY: t_disturb_config
    USE mo_fuel_config_class, ONLY: t_fuel_config

    CHARACTER(len=*), INTENT(in) :: master_namelist_filename

    TYPE(t_sse_config), POINTER :: sse__conf
    TYPE(t_hydro_config), POINTER :: hydro__conf
    TYPE(t_rad_config), POINTER :: rad__conf
    TYPE(t_carbon_config), POINTER :: carbon__conf
    TYPE(t_disturb_config), POINTER :: disturb__conf
    TYPE(t_fuel_config), POINTER :: fuel__conf

    TYPE(t_jsb_model), POINTER :: model => NULL()

    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid

    INTEGER :: no_of_models, model_idx, iproc

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:jsbach_setup_models'

    CALL message(TRIM(routine), 'starting basic setup')

    CALL init_parallel_iface

    CALL init_jsb_master_control(master_namelist_filename)

    !IF (jsbach_runs_standalone()) CALL init_mpi_communicators !TODO-check-with-RMS: remove: done in the driver!

    no_of_models = get_no_of_models()
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models = no_of_models

    ! Create models according to master namelist
    DO model_idx=1,no_of_models
      NULLIFY(model)
      model => new_model(jsb_models_nml(model_idx)%model_id,                &
        &                  jsb_models_nml(model_idx)%model_name,              &
        &                  jsb_models_nml(model_idx)%model_shortname,         &
        &                  jsb_models_nml(model_idx)%model_description,       &
        &                  jsb_models_nml(model_idx)%model_namelist_filename)
      jsbach%models(jsb_models_nml(model_idx)%model_id)%m => model
    END DO

    ! Configure models according to model namelist(s)
    DO model_idx=1,no_of_models
       NULLIFY(model)
       model => jsbach%models(model_idx)%m
       model%config => new_model_config(model%namelist_filename)

       CALL message(TRIM(routine), 'New jsbach model: '//TRIM(model%shortname))
       CALL message('     ', '... ID: '//TRIM(int2string(model%id)))
       CALL message('     ', '... Usecase: '//TRIM(model%config%usecase))
       CALL message('     ', '... Using configuration from namelist: '//TRIM(model%namelist_filename))
       CALL message('     ', '... Using tile fractions from: '//TRIM(model%config%fract_filename))

       IF (jsbach_runs_standalone()) THEN
         CALL init_time(TRIM(model%namelist_filename), TRIM(model%shortname), jsbach_is_restarted())
       ELSE
         CALL init_time()
       END IF

      CALL model%Set_mode(model%config%hsm_mode)

      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(surface)

      ALLOCATE(model%processes(max_no_of_processes))
      DO iproc=1,max_no_of_processes
        model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename))
      END DO

      CALL model%Configure_processes()

      SELECT TYPE (cc => model%processes(sse_)%p%config); TYPE IS (t_sse_config) ; sse__conf => cc ; END SELECT
      SELECT TYPE (cc => model%processes(hydro_)%p%config); TYPE IS (t_hydro_config) ; hydro__conf => cc ; END SELECT
      SELECT TYPE (cc => model%processes(rad_)%p%config); TYPE IS (t_rad_config) ; rad__conf => cc ; END SELECT
      SELECT TYPE (cc => model%processes(carbon_)%p%config); TYPE IS (t_carbon_config) ; carbon__conf => cc ; END SELECT
      SELECT TYPE (cc => model%processes(fuel_)%p%config); TYPE IS (t_fuel_config) ; fuel__conf => cc ; END SELECT
      SELECT TYPE (cc => model%processes(disturb_)%p%config); TYPE IS (t_disturb_config) ; disturb__conf => cc ; END SELECT

      ! Crosscheck process configs
      IF (TRIM(model%config%tpe_scheme) /= '') THEN
        CALL message(TRIM(routine), 'Running terraplanet experiment with scheme '//TRIM(model%config%tpe_scheme))
        IF (model%config%use_lakes) THEN
          CALL message('    ', '... using lakes')
        ELSE
          CALL message('    ', '... not using lakes')
        END IF
        IF (model%config%use_glacier) THEN
          CALL message('    ', '... using glaciers')
        ELSE
          CALL message('    ', '... not using glaciers')
        END IF
      END IF

      IF (model%config%l_compat401) THEN
        CALL message(TRIM(routine), 'Using settings compatible to jsbach4.01:')
        CALL message('    ', '... using soil heat capacity from mineral soil map')
        CALL message('    ', '... using soil heat conductivity from mineral soil map')
        CALL message('    ', '... using fixed snow heat capacity and snow heat conductivity')
        CALL message('    ', '... not using multi-layer snow model')
        CALL message('    ', '... not using dynamical snow density and thickness')
        CALL message('    ', '... not using frozen ground')
        CALL message('    ', '... not using organic soil fractions')
        CALL message('    ', '... using canopy albedo from ini file')
        sse__conf%l_heat_cap_dyn = .FALSE.
        sse__conf%l_heat_cond_dyn = .FALSE.
        sse__conf%l_heat_cap_map = .FALSE.
        sse__conf%l_heat_cond_map = .FALSE.
        sse__conf%l_snow = .FALSE.
        sse__conf%l_dynsnow = .FALSE.
        sse__conf%l_freeze = .FALSE.
        hydro__conf%l_organic = .FALSE.
        rad__conf%use_alb_canopy = .TRUE.
      END IF

      IF (rad__conf%use_alb_canopy .OR. TRIM(model%config%usecase) == 'jsbach_lite') THEN
        CALL message(TRIM(routine), 'Using canopy albedo from ini file '//TRIM(rad__conf%bc_filename)//&
          &                         ' ... therefore no influence of vegetation on albedo')
      ELSE
        CALL message(TRIM(routine), 'Using canopy albedo from lctlib file')
      END IF

      IF (disturb__conf%active .AND. .NOT. fuel__conf%active) THEN
        fuel__conf%active = .TRUE.
        CALL message(TRIM(routine), 'Switching on fuel process since disturbance process is active')
      END IF
      IF (disturb__conf%active .AND. .NOT. carbon__conf%active) THEN
        carbon__conf%active = .TRUE.
        CALL message(TRIM(routine), 'Switching on carbon process since disturbance process is active')
      END IF

    END DO

    IF (debug_on()) CALL message(TRIM(routine), 'done basic setup')

  END SUBROUTINE jsbach_setup_models

  SUBROUTINE jsbach_setup_tiles(model_id)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_model_usecases, ONLY: init_usecase
    USE mo_jsb_lctlib_class,   ONLY: Read_lctlib

    INTEGER, INTENT(in) :: model_id

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:jsbach_setup_tiles'

    model => Get_model(model_id)

    model%lctlib => Read_lctlib('lctlib_nlct21.def')
    CALL message(TRIM(routine), 'Imported lctlib from '//'lctlib_nlct21.def')

    IF (model%config%usecase /= '') THEN
      ! Shortcut model configuration through usecase
      CALL init_usecase(model)
    ELSE
      ! Explicit model configuration through namelist
      ! TBD
    END IF

  END SUBROUTINE jsbach_setup_tiles

#endif
END MODULE mo_jsb_base
