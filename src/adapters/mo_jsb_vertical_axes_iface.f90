!>
!! @brief Contains the adapter for mo_jsb_vertical_axes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2017-03-14)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
MODULE mo_jsb_vertical_axes_iface
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp, dp
  USE mo_exception, ONLY: finish, message

#ifdef __ICON__
  USE mo_name_list_output_zaxes_types, ONLY: t_verticalAxis, t_verticalAxisList
  USE mo_cdi_ids,                      ONLY: set_vertical_grid, t_Vgrid
  USE mo_zaxis_type,                   ONLY: zaxisTypeList
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_verticalAxisList, Setup_jsb_vertical_axis, t_Vgrid, Set_jsb_restart_vgrid

#ifndef __ICON__
  ! Dummy
  TYPE t_verticalAxis
  END TYPE t_verticalAxis
  TYPE t_verticalAxisList
  END TYPE t_verticalAxisList
  TYPE t_Vgrid
  END TYPE
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_vertical_axes_iface'

CONTAINS

#ifdef __ICON__

  SUBROUTINE Setup_jsb_vertical_axis(verticalAxisList, zaxis_id, name, length, longname, units, levels, lbounds, ubounds)

    TYPE(t_verticalAxisList), INTENT(inout) :: verticalAxisList
    INTEGER,                  INTENT(in)    :: zaxis_id
    CHARACTER(len=*),         INTENT(in)    :: name
    INTEGER,                  INTENT(in)    :: length
    CHARACTER(len=*),         INTENT(in)    :: longname
    CHARACTER(len=*),         INTENT(in)    :: units
    REAL(wp),                 INTENT(in)    :: levels(:), lbounds(:), ubounds(:)

    TYPE(t_verticalAxis), POINTER :: axis

    CALL verticalAxisList%append( &
      & t_verticalAxis(zaxisTypeList%getEntry(zaxis_id), &
      &                length,                           &
      &                zaxisLevels   = levels,           &
      &                zaxisLbounds  = lbounds,          &
      &                zaxisUbounds  = ubounds,          &
      &                zaxisUnits    = units,            &
      &                zaxisName     = name,             &
      &                zaxisLongname = longname          &
      &               ) &
      & )
    axis => verticalAxisList%getEntry(icon_zaxis_type=zaxis_id)

  END SUBROUTINE Setup_jsb_vertical_axis

  SUBROUTINE Set_jsb_restart_vgrid(vgrid_defs, count, zaxis_id, levels)

    TYPE(t_Vgrid), INTENT(inout) :: vgrid_defs(:)
    INTEGER,       INTENT(inout) :: count
    INTEGER, VALUE               :: zaxis_id
    REAL(wp),      INTENT(in)    :: levels(:)

    CALL set_vertical_grid(vgrid_defs, count, zaxis_id, levels)

  END SUBROUTINE Set_jsb_restart_vgrid

#else

  SUBROUTINE Setup_jsb_vertical_axis(verticalAxisList, zaxis_id, name, length, longname, units, levels, lbounds, ubounds)
    ! Dummy
    TYPE(t_verticalAxisList), INTENT(inout) :: verticalAxisList
    INTEGER,                  INTENT(in)    :: zaxis_id
    CHARACTER(len=*),         INTENT(in)    :: name
    INTEGER,                  INTENT(in)    :: length
    CHARACTER(len=*),         INTENT(in)    :: longname
    CHARACTER(len=*),         INTENT(in)    :: units
    REAL(wp),                 INTENT(in)    :: levels(:), lbounds(:), ubounds(:)

  END SUBROUTINE Setup_jsb_vertical_axis

  SUBROUTINE Set_jsb_restart_vgrid(vgrid_defs, count, zaxis_id, levels)
    ! Dummy
    TYPE(t_Vgrid), INTENT(inout) :: vgrid_defs(:)
    INTEGER,       INTENT(inout) :: count
    INTEGER, VALUE               :: zaxis_id
    REAL(wp),      INTENT(in)    :: levels(:)
  END SUBROUTINE Set_jsb_restart_vgrid

#endif

#endif
END MODULE mo_jsb_vertical_axes_iface