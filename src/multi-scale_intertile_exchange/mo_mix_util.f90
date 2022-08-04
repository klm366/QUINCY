!>
!! @brief Contains utilities for the the multi-scale intertile exchange process
!!
!! @author
!! Tobias Stacke, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Tobias Stacke (2022-08-03)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_mix_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,  ONLY: ON_TILE_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: find_first_descendant

  INTERFACE find_first_descendant
    MODULE PROCEDURE find_first_descendant
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_util'

CONTAINS

  RECURSIVE SUBROUTINE find_first_descendant(iproc, parenttile, targettile)

    INTEGER, INTENT(in)                    :: iproc
    CLASS(t_jsb_tile_abstract), INTENT(in) :: parenttile !< Parent tile
    CLASS(t_jsb_tile_abstract), POINTER    :: childtile  !< First child 
    CLASS(t_jsb_tile_abstract), POINTER    :: targettile !< First descendant with iproc processes actively computing something

    INTEGER :: n, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':find_first_descendant'

    ! Get first child
    childtile => parenttile%Get_first_child_tile()

    IF (childtile%process_action(iproc) /= ON_TILE_) THEN
      ! No iproc processes running on this tile, search for children
      CALL find_first_descendant(iproc, childtile, targettile)
    ELSE
      targettile => childtile
    END IF

  END SUBROUTINE find_first_descendant

#endif
END MODULE mo_mix_util
