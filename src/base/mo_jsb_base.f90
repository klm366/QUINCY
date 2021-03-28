!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief Contains main JSBACH structure and methods.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!!  First version by Reiner Schnur (2013-04-04)
!!  Modification to fit QUINCY by Jan Engel MPI-BGC (2019-01-08)
!!
MODULE mo_jsb_base
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

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
  !! create the surface vgrid \n
  !! run process configurations, which can create other vgrids \n
  !! JSBACH4: this routine is named jsbach_setup_models()
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE setup_model_vgrid_procconf(namelist_filename)

    USE mo_jsb_process_factory,   ONLY: max_no_of_processes, Create_process
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_config_class,      ONLY: new_model_config
    USE mo_jsb_control,           ONLY: get_no_of_models
    USE mo_jsb_model_class,       ONLY: t_jsb_model, new_model
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,              ONLY: Register_vgrid
    USE mo_jsb_grid_class,        ONLY: ZAXIS_SURFACE   ! in JSBACH4 used from mo_jsb_io    *slmtemp
    ! process config
    dsl4jsb_Use_config(SB_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*), INTENT(in)  :: namelist_filename
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER    :: model => NULL()
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_surface                    ! Vertical grid
    INTEGER                       :: no_of_models, model_idx, iproc
    CHARACTER(len=*), PARAMETER   :: routine = 'mo_jsb_base:setup_model_vgrid_procconf'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    ! note: the according Get_config is called below after the proc_config was created


    !< 0.9
    !! 
    CALL message(TRIM(routine), 'starting basic setup')


    !> 1.0 model
    !!
    !! create just one model and its config
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
    model%config               => new_model_config(namelist_filename)   ! read BASE_CTL from namelist


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
      model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename))
    END DO
    ! run configuration for each process  (mo_jsb_process_class:Configure_process() -> t_PROC_config:Init -> mo_PROC_config_class:Init_PROC_config)
    ! incl. creating vertical grids:
    !   canopy_layer   in ASSIMI_  (ncanopy is part of the &ASSIMILATION_CTL namelist)
    !   soil_layer_sb  in SSE_     (number of soil layers is part of the &SSE_CTL namelist)
    !   spectral_bands in RAD_     (number of spectral bands the light is divided into, vector length is currently hardcoded in rad_config)
    CALL model%Configure_processes()

    !> 3.1 crosscheck process configs
    !!
    dsl4jsb_Get_config(SB_)
    SELECT CASE (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
    CASE ("jsm")
      IF(TRIM(model%config%slm_model_name) /= "soil") THEN
        WRITE(message_text,'(a)') 'setting sb_model_scheme to JSM currently works only with the SOIL standalone model'
        CALL finish(routine, message_text) 
      ENDIF 
    ENDSELECT

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
    USE mo_quincy_set_parameters,     ONLY: parameters_slm_from_file
    USE mo_jsb_math_constants,        ONLY: init_jsb_math_constants
    USE mo_jsb_phys_constants,        ONLY: init_jsb_phys_constants
    USE mo_atmland_constants,         ONLY: init_atmland_constants
    USE mo_rad_constants,             ONLY: init_rad_constants
    USE mo_assimi_constants,          ONLY: init_assimi_constants
    USE mo_pheno_constants,           ONLY: init_pheno_constants
    USE mo_veg_constants,             ONLY: init_veg_constants
    USE mo_sse_constants,             ONLY: init_sse_constants
    USE mo_sb_constants,              ONLY: init_sb_constants

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
    !!
    CALL init_jsb_math_constants      ! may be 1st of the 'init constants'
    CALL init_jsb_phys_constants      ! may be 2nd of the 'init constants'
    CALL init_rad_constants           ! may be 3rd of the 'init constants'
    CALL init_atmland_constants       ! uses hourday                  defined in mo_jsb_math_constants
    CALL init_assimi_constants        ! uses molar_mass_N, Dwv, Dco2  defined in mo_jsb_phys_constants
    CALL init_veg_constants           ! uses Tzero                    defined in mo_jsb_phys_constants
    CALL init_sse_constants
    CALL init_sb_constants
    CALL init_pheno_constants         ! uses Tzero                    defined in mo_jsb_phys_constants


    !> 2.0 parameter modification for parameter sensitivity analysis
    !! 
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
    USE mo_quincy_set_parameters,     ONLY: set_parameter_lctlib_values_with_list_values

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
    model%lctlib => Read_lctlib('lctlib_slm_nlct8.def')
    CALL message(TRIM(routine), 'Imported lctlib from '//'lctlib_slm_nlct8.def')


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
