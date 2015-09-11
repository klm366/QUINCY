!>
!! @brief Contains process class
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2015-09-01)
!!
MODULE mo_jsb_process_class
#ifndef __NO_JSBACH__

  USE mo_jsb_config_class, ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_process

  TYPE t_jsb_process
    INTEGER                      :: id = 0
    CLASS(t_jsb_config), POINTER :: config => NULL()
  CONTAINS
    PROCEDURE :: Configure => Configure_process
  END TYPE t_jsb_process

CONTAINS

  SUBROUTINE Configure_process(this, iproc, namelist_filename)

    USE mo_jsb_process_factory, ONLY: Create_process_config

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: iproc
    CHARACTER(LEN=*),     INTENT(in)    :: namelist_filename

    this%id       = iproc
    this%config  => Create_process_config(iproc, namelist_filename)

  END SUBROUTINE Configure_process

#endif
END MODULE mo_jsb_process_class
