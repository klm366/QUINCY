!> Contains main JSBACH structure and methods.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_base
#ifndef __NO_JSBACH__

!!$  USE mo_jsb_model_class, ONLY: t_jsb_model_m
!!$  USE mo_jsb_class, ONLY: t_jsb, jsbach
  USE mo_jsb_class,   ONLY: jsbach          ! USE here to be able to define as public
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message, message_text, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach
  PUBLIC :: setup_model_vgrid_procconf, quincy_proc_parameter_init, jsbach_setup_tiles

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> set up the model with vgrid and process configs
  !!
  !! create/allocate models \n
  !!  plus the bgcm_store \n
  !! create the surface vgrid \n
  !! run process configurations, which can create other vgrids \n
  !! JSBACH4: this routine is named jsbach_setup_models()
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE setup_model_vgrid_procconf(namelist_filename)

    USE mo_jsb_process_factory,   ONLY: max_no_of_processes, Create_process
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SPQ_, HYDRO_, VEG_, HD_, RAD_, Q_ASSIMI_, Q_PHENO_, SB_, PHYD_
    USE mo_jsb_config_class,      ONLY: new_model_config
    USE mo_jsb_control,           ONLY: get_no_of_models
    USE mo_jsb_model_class,       ONLY: t_jsb_model, new_model
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,              ONLY: Register_vgrid
    USE mo_jsb_io,                ONLY: ZAXIS_SURFACE
    USE mo_lnd_bgcm_idx,          ONLY: set_bgcm_element_idx

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*), INTENT(in)  :: namelist_filename
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER    :: model => NULL()
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_surface                    ! Vertical grid
    INTEGER                       :: no_of_models, model_idx, iproc
    INTEGER                       :: model_scheme                     !< MODEL_JSBACH, MODEL_QUINCY
    CHARACTER(len=*), PARAMETER   :: routine = 'mo_jsb_base:setup_model_vgrid_procconf'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    ! note: the according Get_config is called below after the proc_config was created


    !>0.9
    !>
    CALL message(TRIM(routine), 'starting basic setup')


    !>1.0 model
    !>
    !! create just one model and its config
    !! plus the bgcm_store and related stuff
    !!
    no_of_models        =  get_no_of_models()         ! hardcoded to be 1
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models =  no_of_models
    model_idx           =  1
    model => new_model(model_idx, &                   ! model_id, name, shortname, description, namelist_filename
                       "quincy model", &
                       "quincy", &
                       "Quantifying interactions between terrestrial nutrient cycles and the climate system", &
                       namelist_filename)
    jsbach%models(model_idx)%m => model
    model%config               => new_model_config(namelist_filename)   !< read BASE_CTL from namelist
    model_scheme = model%config%model_scheme                            !< MODEL_QUINCY = 3
    CALL set_bgcm_element_idx(model%config%elements_index_map)          !< set the bgcm index var used in the scientific routines

    !> 2.0 surface vgrid
    !!
    !! create the only vgrid that is definitely needed by the model
    vgrid_surface => new_vgrid('surface', ZAXIS_SURFACE, 1)  ! vertical grid with one layer - needed for all 2D variables
    CALL register_vgrid(vgrid_surface)


    !> 3.0 process configs and other vgrids
    !!
    !! allocate an object for each process \n
    !! run process configurations which might include creating more vgrids
    ALLOCATE(model%processes(max_no_of_processes))
    DO iproc=1,max_no_of_processes
      model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename), model_scheme)
    END DO

    ! point from each associated process%config to the model_config
    DO iproc = 1,max_no_of_processes
      IF (ASSOCIATED(model%processes(iproc)%p)) THEN
        IF (ASSOCIATED(model%processes(iproc)%p%config)) THEN
          model%processes(iproc)%p%config%model_config => model%config
        ENDIF
      ENDIF
    ENDDO

    ! run configuration for each process  (mo_jsb_process_class:Configure_process() -> t_PROC_config:Init -> mo_PROC_config_class:Init_PROC_config)
    ! incl. creating vertical grids:
    !   q_canopy_layer  in Q_ASSIMI_ (ncanopy is part of the &lnd_q_assimi_nml namelist)
    !   quincy_nsoil_w  in HYDRO_    (number of soil layers for water calculations is part of the &jsb_hydro_nml namelist)
    !   quincy_nsoil_e  in SSE_      (number of soil layers for nergy calculations is part of the &jsb_sse_nml namelist)
    !   spectral_bands  in RAD_      (number of spectral bands the light is divided into, vector length is currently hardcoded in rad_config)
    CALL model%Configure_processes()

    IF (debug_on()) CALL message(TRIM(routine), 'done basic setup')

  END SUBROUTINE setup_model_vgrid_procconf


  !-----------------------------------------------------------------------------------------------------
  !> init processes parameters
  !!
  !! must be done prior to jsbach_setup_tiles() which runs Read_lctlib() \n
  !! may have no dependencies to any other than parameters from other processes
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE quincy_proc_parameter_init(model_id)

    USE mo_jsb_model_class,           ONLY: t_jsb_model
    USE mo_jsb_class,                 ONLY: Get_model
    USE mo_qs_set_parameters,         ONLY: parameters_slm_from_file
    USE mo_atmland_constants,         ONLY: init_atmland_constants
    USE mo_qs_atmland_constants,      ONLY: init_qs_atmland_constants
    USE mo_rad_parameters,            ONLY: init_rad_parameters
    USE mo_q_assimi_parameters,       ONLY: init_q_assimi_parameters
    USE mo_q_pheno_parameters,        ONLY: init_q_pheno_parameters
    USE mo_phyd_parameters,           ONLY: init_phyd_parameters
    USE mo_dfire_constants,           ONLY: init_dfire_constants
    USE mo_veg_constants,             ONLY: init_veg_constants
    USE mo_spq_constants,             ONLY: init_spq_constants
    USE mo_sb_constants,              ONLY: init_sb_constants
    USE mo_phyd_constants,            ONLY: init_phyd_constants
    USE mo_q_agr_constants,           ONLY: init_q_agr_constants
    USE mo_turb_constants,            ONLY: init_turb_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER, INTENT(in) :: model_id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:quincy_proc_parameter_init'
    ! ---------------------------
    ! 0.5 Get Memory
    model => Get_model(model_id)



    !> 1.0 parameter init
    !>
    CALL init_rad_parameters
    CALL init_atmland_constants         ! uses hourday                  defined in mo_jsb_math_constants
    CALL init_qs_atmland_constants
    CALL init_q_assimi_parameters
    CALL init_veg_constants             ! uses Tzero                    defined in mo_jsb_physical_constants
    CALL init_spq_constants
    CALL init_sb_constants
    CALL init_q_pheno_parameters
    CALL init_phyd_parameters
    CALL init_dfire_constants
    CALL init_phyd_constants
    CALL init_q_agr_constants
    CALL init_turb_constants


    !> 2.0 parameter modification for parameter sensitivity analysis
    !>
    IF (model%config%set_parameter_values_from_file) THEN
      CALL parameters_slm_from_file%set_param(model)
    ENDIF

  END SUBROUTINE quincy_proc_parameter_init


  !-----------------------------------------------------------------------------------------------------
  !> set up tiles
  !!
  !! read lctlib and modify lctlib-parameter values (if enabled) \n
  !! init usecase
  !!
  !! note: quincy_proc_parameter_init() must run prior to this routine
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE jsbach_setup_tiles(model_id)

    USE mo_jsb_model_class,           ONLY: t_jsb_model
    USE mo_jsb_class,                 ONLY: Get_model
    USE mo_jsb_model_usecases,        ONLY: init_usecase
    USE mo_jsb_lctlib_class,          ONLY: Read_lctlib
    USE mo_qs_set_parameters,     ONLY: set_parameter_lctlib_values_with_list_values

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER, INTENT(in) :: model_id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:jsbach_setup_tiles'
    ! ---------------------------
    ! 0.5 Get Memory
    model => Get_model(model_id)



    !< 1.0 read lctlib
    !!
    model%lctlib => Read_lctlib('lctlib_quincy_nlct14.def')
    CALL message(TRIM(routine), 'Imported lctlib from '//'lctlib_quincy_nlct14.def')


    !> 2.0 lcltib-parameter modification for parameter sensitivity analysis
    !!
    IF (model%config%set_parameter_values_from_file) CALL set_parameter_lctlib_values_with_list_values(model)


    !< 3.0 init usecase
    !!
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
