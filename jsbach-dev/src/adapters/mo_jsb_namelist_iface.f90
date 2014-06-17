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
#ifndef __NO_JSBACH__

  USE mo_exception, ONLY: finish

#ifdef __ICON__
  USE mo_namelist, ONLY: open_nml_icon => open_nml, close_nml_icon => close_nml, &
                         position_nml_icon => position_nml 
  USE mo_io_units, ONLY: nnml
#else
  USE mo_namelist, ONLY: open_nml, close_nml, position_nml_echam => position_nml
#endif
  USE mo_namelist, ONLY: POSITIONED

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
    open_nml = 1

  END FUNCTION open_nml

  !>
  !! Wrapper function to make ICON subroutine close_nml a function take unit as argument,
  !! as in ECHAM. Should be replaced when close_nml in ICON gets updated.
  SUBROUTINE close_nml(nml_handler)

    INTEGER, INTENT(in) :: nml_handler

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':close_nml'

    IF (nml_handler /= 1) CALL finish(TRIM(routine), 'Cannot close namelist unit')

    CALL close_nml_icon()

  END SUBROUTINE close_nml
#endif

  !>
  !! Wrapper function to make ICON subroutine position_nml a function, as in ECHAM.
  !! Should be replaced when open_nml in ICON gets updated.
  INTEGER FUNCTION position_nml(name, nml_handler, lrewind, status)

    CHARACTER(len=*), INTENT(in)            :: name         ! namelist group name
    INTEGER,          INTENT(in)            :: nml_handler  ! file handler
    LOGICAL,          INTENT(in)  ,OPTIONAL :: lrewind      ! default: true
    INTEGER,          INTENT(out) ,OPTIONAL :: status       ! error return value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':position_nml'

#ifdef __ICON__
    CALL position_nml_icon(name, lrewind=lrewind, status=status)

    position_nml = nnml
#else
    position_nml = position_nml_echam(name, nml_handler, rewind=lrewind, status=status)
#endif

  END FUNCTION position_nml

#endif
END MODULE mo_jsb_namelist_iface
