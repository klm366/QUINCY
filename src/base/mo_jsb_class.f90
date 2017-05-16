!>
!! @brief basic container (structure) for models
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg (Jan Engel MPI-BGC thinks Reiner wrote this module)
!!
!! @par Revision History
!!  Reiner Schnur, MPI-M Hamburg
!!  Jan Engel, MPI-BGC Jena: initial integration into SLM and adaption to SLM (20161025) & update 20170516
!!
MODULE mo_jsb_class
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class, ONLY: t_jsb_model, t_jsb_model_m

  IMPLICIT NONE
  PRIVATE

  ! t_jsb needs to be public for the NEC compiler
  PUBLIC :: t_jsb, jsbach, get_model

  TYPE t_jsb
     INTEGER :: no_of_models
     TYPE(t_jsb_model_m), POINTER :: models(:)
  END TYPE t_jsb

  TYPE(t_jsb) :: jsbach

  INTERFACE get_model
     MODULE PROCEDURE get_model_by_id
     !MODULE PROCEDURE get_model_by_name
  END INTERFACE get_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_class'

CONTAINS

!!$  FUNCTION get_grid_by_id(id) RESULT(grid)
!!$
!!$    USE mo_jsb_grid_class, ONLY: t_jsb_grid
!!$
!!$    INTEGER, INTENT(in)        :: id
!!$    TYPE(t_jsb_grid), POINTER :: grid
!!$
!!$    grid => jsbach%models(id)%m%grid
!!$
!!$  END FUNCTION get_grid_by_id

  FUNCTION get_model_by_id(id) RESULT(model)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)        :: id
    TYPE(t_jsb_model), POINTER :: model
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_model_by_id'

    model => jsbach%models(id)%m

  END FUNCTION get_model_by_id

#endif
END MODULE mo_jsb_class
