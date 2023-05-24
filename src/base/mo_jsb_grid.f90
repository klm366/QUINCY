!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_grid
#ifndef __NO_JSBACH__

  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_jsb_grid_class, ONLY: t_jsb_grid, t_jsb_grid_ptr

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_grid, Get_grid

  TYPE(t_jsb_grid_ptr),  ALLOCATABLE, SAVE :: jsbach_grids(:)

CONTAINS

  SUBROUTINE Register_grid(grid)

    TYPE(t_jsb_grid),  POINTER, INTENT(in)    :: grid

    INTEGER :: i, n
    TYPE(t_jsb_grid_ptr), ALLOCATABLE :: grids_temp(:)

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
      n = SIZE(jsbach_grids)
      ALLOCATE(grids_temp(n+1))
      grids_temp(1:n) = jsbach_grids
      grids_temp(n+1) = t_jsb_grid_ptr(grid)
      CALL MOVE_ALLOC(grids_temp, jsbach_grids)
    ELSE
      ALLOCATE(jsbach_grids(1))
      jsbach_grids(1) = t_jsb_grid_ptr(grid)
    END IF

    CALL message(TRIM(routine), 'registered new JSBACH grid instance '//TRIM(grid%name))

    !$ACC ENTER DATA PCREATE(grid, grid%lat, grid%lon, grid%area, grid%lsm, grid%lsf)
    !$ACC UPDATE DEVICE(grid%lat, grid%lon, grid%area, grid%lsm, grid%lsf)

  END SUBROUTINE Register_grid

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

#endif
END MODULE mo_jsb_grid
