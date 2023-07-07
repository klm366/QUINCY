!>
!! @brief Contains methods for JSBACH finalization.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                          by Reiner Schnur (2018-06-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_model_final
#ifndef __NO_JSBACH__

  ! USE mo_kind,             ONLY: wp
  ! USE mo_exception,        ONLY: message, finish

  ! USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_control,      ONLY: l_timer_host, timer_on !,get_no_of_models, debug_on
  ! USE mo_jsb_model_class,  ONLY: t_jsb_model

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_finalize

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_final'

CONTAINS

  SUBROUTINE jsbach_finalize

    USE mo_timer, ONLY: print_timer

    ! TYPE(t_jsb_model), POINTER  :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_finalize'

    IF (.NOT. l_timer_host .AND. timer_on()) CALL print_timer()

  END SUBROUTINE jsbach_finalize

#endif
END MODULE mo_jsb_model_final
