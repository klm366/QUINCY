MODULE mo_jsb_grid
#ifndef __NO_JSBACH__

  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_jsb_grid_class, ONLY: t_jsb_grid, t_jsb_grid_ptr, t_jsb_vgrid, t_jsb_vgrid_ptr

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_grid, Register_vgrid, Get_grid, Get_vgrid

  TYPE(t_jsb_grid_ptr),  ALLOCATABLE, SAVE :: jsbach_grids(:)

CONTAINS

  SUBROUTINE Register_grid(grid)

    USE mo_jsb_vertical_axes, ONLY: jsbach_vgrids

    TYPE(t_jsb_grid),  POINTER, INTENT(in)    :: grid

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Register_grid'

    IF (ALLOCATED(jsbach_grids)) THEN
      DO i=1,SIZE(jsbach_grids)
        IF (ASSOCIATED(jsbach_grids(i)%p)) THEN
          IF (jsbach_grids(i)%p%id == grid%id) THEN
            RETURN ! grid already exists and is already registered
          END IF
        END IF
      END DO
      ! grid is new
      jsbach_grids = (/jsbach_grids, t_jsb_grid_ptr(grid)/)
    ELSE
      jsbach_grids = (/t_jsb_grid_ptr(grid)/)
    END IF

    CALL message(TRIM(routine), 'registered new JSBACH grid instance '//TRIM(grid%name))

  END SUBROUTINE Register_grid

  SUBROUTINE Register_vgrid(vgrid)

    USE mo_jsb_io_iface, ONLY: Create_zaxis
    USE mo_jsb_vertical_axes, ONLY: jsbach_vgrids

    TYPE(t_jsb_vgrid), POINTER, INTENT(inout) :: vgrid

    INTEGER :: i, length

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Register_vgrid'

    IF (ALLOCATED(jsbach_vgrids)) THEN
      DO i=1,SIZE(jsbach_vgrids)
        IF (ASSOCIATED(jsbach_vgrids(i)%p)) THEN
          IF (jsbach_vgrids(i)%p%ZaxisID == vgrid%ZaxisID) THEN
            RETURN ! axis already exists and is already registered
          END IF
        END IF
      END DO
      ! Axis is new
      jsbach_vgrids = (/jsbach_vgrids, t_jsb_vgrid_ptr(vgrid)/)
    ELSE
      jsbach_vgrids = (/t_jsb_vgrid_ptr(vgrid)/)
    END IF
    vgrid%id = SIZE(jsbach_vgrids)

    length = vgrid%n_levels
    CALL Create_zaxis(vgrid%ZaxisID, vgrid%echamZaxisIdx, &
      & vgrid%cdi_axis_type, length, vgrid%name, vgrid%levels(1:length), vgrid%units, vgrid%longname, &
      & vgrid%lbounds(1:length), vgrid%ubounds(1:length))

    WRITE(message_text, '(A,A,2I3)') 'Registered new vgrid: name, ID, length= ',TRIM(vgrid%name), vgrid%ZaxisID, length
    CALL message(TRIM(routine), message_text)

  END SUBROUTINE Register_vgrid

  FUNCTION Get_grid(id) RESULT(grid)

    INTEGER, INTENT(in)        :: id
    TYPE(t_jsb_grid), POINTER  :: grid

    INTEGER :: i
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Get_grid'

    grid => NULL()
    DO i=1,SIZE(jsbach_grids)
      IF (ASSOCIATED(jsbach_grids(i)%p)) THEN
        IF (jsbach_grids(i)%p%id == id) THEN
          grid => jsbach_grids(i)%p
          EXIT
        END IF
      END IF
    END DO

    IF (.NOT. ASSOCIATED(grid)) THEN
      CALL finish(TRIM(routine), 'Grid not found.')
    END IF

  END FUNCTION Get_grid

  FUNCTION Get_vgrid(name) RESULT(vgrid)

    USE mo_jsb_vertical_axes, ONLY: jsbach_vgrids

    CHARACTER(len=*),           INTENT(in)    :: name
    TYPE(t_jsb_vgrid), POINTER                :: vgrid

    INTEGER :: i
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Get_vgrid'

    vgrid => NULL()
    DO i=1,SIZE(jsbach_vgrids)
      IF (ASSOCIATED(jsbach_vgrids(i)%p)) THEN
        IF (jsbach_vgrids(i)%p%id > 0 .AND. TRIM(jsbach_vgrids(i)%p%name) == TRIM(name)) THEN
          vgrid => jsbach_vgrids(i)%p
          EXIT
        END IF
      END IF
    END DO

    IF (.NOT. ASSOCIATED(vgrid)) THEN
      CALL finish(TRIM(routine), 'Vgrid '//TRIM(name)//' not found.')
    END IF

  END FUNCTION Get_vgrid

#endif
END MODULE mo_jsb_grid
