!>
!! @brief Contains interface to namelist handling
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-02-06)
!!
MODULE mo_jsb_namelist_iface

  USE mo_exception, ONLY: finish

#ifdef __ICON__
  USE mo_namelist, ONLY: open_nml_icon => open_nml, close_nml_icon => close_nml
  USE mo_io_units, ONLY: nnml
#else
  USE mo_namelist, ONLY: open_nml, close_nml
#endif
  USE mo_namelist, ONLY: POSITIONED, position_nml

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: open_nml, POSITIONED, position_nml, close_nml

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_namelist_iface'

CONTAINS

#ifdef __ICON__
  !>
  !! Wrapper function to make ICON subroutine open_nml a function, as in ECHAM.
  !! Should be replaced when open_nml in ICON gets updated.
  INTEGER FUNCTION open_nml(namelist_filename)

    CHARACTER(len=*), INTENT(in) :: namelist_filename

    LOGICAL :: unit_is_occupied

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_nml'

    INQUIRE(UNIT=nnml, OPENED=unit_is_occupied)
    IF (unit_is_occupied) CALL finish(TRIM(routine), 'namelist unit is already open')

    CALL open_nml_icon(TRIM(namelist_filename))
    open_nml = nnml

  END FUNCTION open_nml

  !>
  !! Wrapper function to make ICON subroutine close_nml a function take unit as argument,
  !! as in ECHAM. Should be replaced when close_nml in ICON gets updated.
  SUBROUTINE close_nml(nml_unit)

    INTEGER, INTENT(in) :: nml_unit

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':close_nml'

    IF (nml_unit /= nnml) CALL finish(TRIM(routine), 'Cannot close namelist unit')

    CALL close_nml_icon()

  END SUBROUTINE close_nml
#endif

END MODULE mo_jsb_namelist_iface
