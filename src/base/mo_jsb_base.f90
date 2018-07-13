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
  USE mo_exception,   ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach
  PUBLIC :: init_base !, register_models

CONTAINS

  SUBROUTINE init_base(master_namelist_filename)

    USE mo_jsb_control,         ONLY: get_no_of_models
    USE mo_jsb_model_class,     ONLY: t_jsb_model, new_model

    CHARACTER(len=*), INTENT(in), OPTIONAL :: master_namelist_filename

    TYPE(t_jsb_model), POINTER :: current => NULL()

    INTEGER :: no_of_models, model_id

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:init_base'

    CALL message(TRIM(routine), 'starting basic initialization')

    no_of_models = get_no_of_models()
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models = no_of_models

    model_id = 1
      current => new_model(model_id, "slm_nm", "slm_snm", "slm_descr", "namelist.slm") ! model_id, name, shortname, description, namelist_filename
      jsbach%models(model_id)%m => current

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
