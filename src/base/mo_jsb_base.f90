!>
!! @brief Contains main JSBACH structure and methods.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_base
#ifndef __NO_JSBACH__

!!$  USE mo_jsb_model_class, ONLY: t_jsb_model_m
!!$  USE mo_jsb_class, ONLY: t_jsb, jsbach
  USE mo_jsb_class,   ONLY: jsbach
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message_text, message, finish
  USE mo_util,        ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach
  PUBLIC :: init_base !, register_models

CONTAINS

  SUBROUTINE init_base(master_namelist_filename)

    USE mo_jsb_control,        ONLY: jsb_models_nml, get_no_of_models, init_jsb_master_control, &
      & jsbach_runs_standalone, jsbach_is_restarted
    USE mo_jsb_parallel,       ONLY: init_parallel_iface, init_mpi_communicators
    USE mo_jsb_model_class,    ONLY: t_jsb_model, new_model
    USE mo_jsb_config_class,   ONLY: new_model_config
    USE mo_jsb_time,           ONLY: init_time
    !USE mo_jsb_domain,         ONLY: init_domains

    CHARACTER(len=*), INTENT(in) :: master_namelist_filename

    TYPE(t_jsb_model), POINTER :: current => NULL()

    INTEGER :: no_of_models, model_idx

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:init_base'

    CALL message(TRIM(routine), 'starting basic initialization')

    CALL init_jsb_master_control(master_namelist_filename)

    CALL init_parallel_iface

    !IF (jsbach_runs_standalone()) CALL init_mpi_communicators !TODO-check-with-RMS: remove: done in the driver!

    no_of_models = get_no_of_models()
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models = no_of_models

    ! Create models according to master namelist
    DO model_idx=1,no_of_models
      NULLIFY(current)
      current => new_model(jsb_models_nml(model_idx)%model_id,                &
        &                  jsb_models_nml(model_idx)%model_name,              &
        &                  jsb_models_nml(model_idx)%model_shortname,         &
        &                  jsb_models_nml(model_idx)%model_description,       &
        &                  jsb_models_nml(model_idx)%model_namelist_filename)
      jsbach%models(jsb_models_nml(model_idx)%model_id)%m => current
    END DO

    ! Configure models according to model namelist(s)
    DO model_idx=1,no_of_models
       NULLIFY(current)
       current => jsbach%models(model_idx)%m
       current%config => new_model_config(current%namelist_filename)

       CALL message(TRIM(routine), 'New jsbach model: '//TRIM(current%shortname))
       CALL message('     ', '... ID: '//TRIM(int2string(current%id)))
       CALL message('     ', '... Usecase: '//TRIM(current%config%usecase))
       CALL message('     ', '... Reading configuration from namelist: '//TRIM(current%namelist_filename))
       CALL message('     ', '... Reading tile fractions from: '//TRIM(current%config%fract_filename))

       IF (jsbach_runs_standalone()) THEN
         CALL init_time(TRIM(current%namelist_filename), TRIM(current%shortname), jsbach_is_restarted())
       ELSE
         CALL init_time()
       END IF

    END DO

    !IF (jsbach_runs_standalone) CALL init_domains(jsbach%models(:)%grid%filename) !TODO-check-with-RMS: delete!?

    IF (debug_on()) CALL message(TRIM(routine), 'done basic initialization')

  END SUBROUTINE init_base

!!$  SUBROUTINE register_model(model)
!!$
!!$    USE mo_jsb_model_class, ONLY: t_jsb_model
!!$
!!$#ifdef HAVE_F2003
!!$    TYPE(t_jsb_model), POINTER, INTENT(in) :: model
!!$#else
!!$    TYPE(t_jsb_model), POINTER             :: model
!!$#endif
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:register_model'
!!$
!!$    CALL message(TRIM(routine), 'going to register new JSBACH model instance')
!!$
!!$    !jsbach%nmodels = jsbach%nmodels + 1
!!$    model%id = jsbach%nmodels
!!$    land%models(jsbach%nmodels)%m => model
!!$
!!$    CALL message(TRIM(routine), 'registered new JSBACH model instance '//model%name)
!!$
!!$  END SUBROUTINE register_model

!!$  FUNCTION get_model_by_id(id) RESULT(model)
!!$
!!$    USE mo_jsb_model_class, ONLY: t_jsb_model
!!$
!!$    INTEGER, INTENT(in)        :: id
!!$    TYPE(t_jsb_model), POINTER :: model
!!$
!!$    model => jsbach%models(id)%m
!!$
!!$  END FUNCTION get_model_by_id

#endif
END MODULE mo_jsb_base
