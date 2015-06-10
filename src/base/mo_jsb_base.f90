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
  USE mo_exception,   ONLY: message_text, message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach
  PUBLIC :: init_base !, register_model
  PUBLIC :: get_model, get_tile

!!$  TYPE t_jsb
!!$    INTEGER :: no_of_models
!!$    TYPE(t_jsb_model_m), POINTER :: models(:)
!!$  END TYPE t_jsb


!!$  TYPE(t_jsb) :: jsbach

  INTERFACE get_model
     MODULE PROCEDURE get_model_by_id
     !MODULE PROCEDURE get_model_by_name
  END INTERFACE get_model

  INTERFACE get_tile
    MODULE PROCEDURE get_tile_by_name
  END INTERFACE get_tile

CONTAINS

  SUBROUTINE init_base(master_namelist_filename)

    USE mo_jsb_control,        ONLY: jsb_models_nml, get_no_of_models, init_jsb_master_control, jsbach_runs_standalone
    USE mo_jsb_parallel_iface, ONLY: init_parallel_iface, init_mpi_communicators
    USE mo_jsb_model_class,    ONLY: t_jsb_model, new_model
    USE mo_jsb_config_class,   ONLY: new_model_config
    !USE mo_jsb_domain,         ONLY: init_domains

    CHARACTER(len=*), INTENT(in) :: master_namelist_filename

    TYPE(t_jsb_model), POINTER :: current => NULL()

    INTEGER :: no_of_models, model_idx

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_base:init_base'

    CALL message(TRIM(routine), 'starting basic initialization')

    CALL init_jsb_master_control(master_namelist_filename)

    CALL init_parallel_iface

    IF (jsbach_runs_standalone()) CALL init_mpi_communicators

    no_of_models = get_no_of_models()
    ALLOCATE(jsbach%models(no_of_models))
    jsbach%no_of_models = no_of_models

    ! Create models according to master namelist
    DO model_idx=1,no_of_models
      NULLIFY(current)
      current => new_model(jsb_models_nml(model_idx)%model_name,              &
                           jsb_models_nml(model_idx)%model_shortname,         &
                           jsb_models_nml(model_idx)%model_description,       &
                           jsb_models_nml(model_idx)%model_namelist_filename)
      current%id = model_idx
      jsbach%models(model_idx)%m => current
    END DO

    ! Configure models according to model namelist(s)
    DO model_idx=1,no_of_models
       NULLIFY(current)
       current => jsbach%models(model_idx)%m
       current%config => new_model_config(current%namelist_filename)
    END DO

    !IF (jsbach_runs_standalone) CALL init_domains(jsbach%models(:)%grid%filename)

    CALL message(TRIM(routine), 'done basic initialization')

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

  FUNCTION get_model_by_id(id) RESULT(model)

    USE mo_jsb_model_class, ONLY: t_jsb_model

    INTEGER, INTENT(in)        :: id
    TYPE(t_jsb_model), POINTER :: model

    model => jsbach%models(id)%m

  END FUNCTION get_model_by_id

  RECURSIVE FUNCTION get_tile_by_name(tile, name) RESULT(found)

    USE mo_jsb_tile_class, ONLY: t_jsb_tile

    TYPE(t_jsb_tile), POINTER    :: tile
    CHARACTER(len=*), INTENT(in) :: name
    TYPE(t_jsb_tile), POINTER    :: found, current

    INTEGER :: i

    found => NULL()

    IF (TRIM(tile%name) == TRIM(name)) THEN
      found => tile
    ELSE IF (tile%no_of_sub_tiles > 0) THEN
      DO i=1,tile%no_of_sub_tiles
        current => tile%sub_tiles(i)
        found => get_tile_by_name(current, name)
        IF (ASSOCIATED(found)) EXIT
      END DO
    END IF

  END FUNCTION get_tile_by_name

#endif
END MODULE mo_jsb_base
