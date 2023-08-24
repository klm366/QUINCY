!>
!! @brief Contains JSBACH vertical axes structure, registry and methods
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

MODULE mo_jsb_vertical_axes
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: message, message_text, finish
  USE mo_jsb_vertical_axes_iface, ONLY: t_verticalAxisList, Setup_jsb_vertical_axis, &
    &                                   t_Vgrid, Set_jsb_restart_vgrid

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_vgrid, t_jsb_vgrid_ptr, new_vgrid, Register_vgrid, Get_vgrid, delete_vgrid, Get_no_of_vgrids
  PUBLIC :: setup_zaxes_jsbach, set_vertical_grids_jsbach

  TYPE t_jsb_vgrid
    INTEGER ::           idx          = 0  !< Index of vgrid
    INTEGER ::           id           = 0  !< ID of vgrid (according to ENUMERATOR)
    CHARACTER(len=20) :: name         = ''
    CHARACTER(len=50) :: longname     = ''
    CHARACTER(len=10) :: units        = ''
    INTEGER :: n_levels               = 0
    REAL(wp), ALLOCATABLE :: levels(:)
    REAL(wp), ALLOCATABLE :: lbounds(:)
    REAL(wp), ALLOCATABLE :: ubounds(:)
    REAL(wp), ALLOCATABLE :: dz(:)
    INTEGER :: cdi_axis_type          = 0
    INTEGER :: ZaxisID                = 0   ! ICON: zaxis type (ZA_*), ECHAM: cdi zaxis id
    INTEGER :: echamZaxisIdx          = 0   ! ICON: N/A, ECHAM: levelindx for dimension
  END TYPE t_jsb_vgrid

  TYPE t_jsb_vgrid_ptr
    TYPE(t_jsb_vgrid), POINTER :: p => NULL()
  END TYPE t_jsb_vgrid_ptr

  ENUM, BIND(C)
    ENUMERATOR :: VGRID_SURFACE=1, VGRID_JSB_WATER, VGRID_JSB_ENERGY, VGRID_JSB_SNOW, VGRID_JSB_CANOPY, VGRID_JSB_PFTS, &
      &           VGRID_AGE_CLASSES,                                                                    &
      &           VGRID_HD_OVERLFLOW, VGRID_HD_BASEFLOW, VGRID_HD_RIVERFLOW,                            &
      &           VGRID_Q_SOIL, VGRID_Q_CANOPY,                                                         &
      &           INVALID
  END ENUM
  PUBLIC :: VGRID_SURFACE, VGRID_JSB_WATER, VGRID_JSB_ENERGY, VGRID_JSB_SNOW, VGRID_JSB_CANOPY, VGRID_JSB_PFTS, &
    &       VGRID_AGE_CLASSES,                                                                  &
    &       VGRID_HD_OVERLFLOW, VGRID_HD_BASEFLOW, VGRID_HD_RIVERFLOW,                          &
    &       VGRID_Q_SOIL, VGRID_Q_CANOPY,                                                       &
    &       INVALID

  TYPE(t_jsb_vgrid_ptr), ALLOCATABLE, SAVE :: jsbach_vgrids(:)

  INTERFACE Get_vgrid
    MODULE PROCEDURE Get_vgrid_by_name
    MODULE PROCEDURE Get_vgrid_by_idx_or_id
  END INTERFACE Get_vgrid

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_vertical_axes'

CONTAINS

  FUNCTION new_vgrid(name, id, cdi_axis_type, length, longname, units, levels, lbounds, ubounds, dz) RESULT(new)


    CHARACTER(len=*),           INTENT(in) :: name
    INTEGER,                    INTENT(in) :: id
    INTEGER,                    INTENT(in) :: cdi_axis_type
    INTEGER,                    INTENT(in) :: length
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: longname
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: units
    REAL(wp),         OPTIONAL, INTENT(in) :: levels(length)
    REAL(wp),         OPTIONAL, INTENT(in) :: lbounds(length)
    REAL(wp),         OPTIONAL, INTENT(in) :: ubounds(length)
    REAL(wp),         OPTIONAL, INTENT(in) :: dz(length)

    TYPE(t_jsb_vgrid), POINTER :: new

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_vgrid'

    IF (TRIM(name) == '') CALL finish(routine, 'Must specify axis name')
    IF (id > INVALID) CALL finish(routine, 'Vertical grid ID not known ... register it in '//routine//' first.')
    IF (length < 1) CALL finish(routine, TRIM(name)//' - length < 1')

    ALLOCATE(new)

    new%cdi_axis_type = cdi_axis_type
    new%name = TRIM(name)
    new%id = id
    new%n_levels = length

    IF (PRESENT(longname)) THEN
      new%longname = TRIM(longname)
    ELSE
      new%longname = ''
    END IF
    IF (PRESENT(units)) THEN
      new%units  = TRIM(units)
    ELSE
      new%units = ''
    END IF

    ALLOCATE(new%levels(length), new%lbounds(length), new%ubounds(length), new%dz(length))
    new%levels (:) = 0._wp
    new%lbounds(:) = 0._wp
    new%ubounds(:) = 0._wp
    new%dz     (:) = 0._wp
    !$ACC ENTER DATA PCREATE(new%levels, new%lbounds, new%ubounds, new%dz)

    IF (PRESENT(levels))  new%levels (:) = levels(:)
    IF (PRESENT(lbounds)) new%lbounds(:) = lbounds(:)
    IF (PRESENT(ubounds)) new%ubounds(:) = ubounds(:)
    IF (PRESENT(dz)) THEN
      IF (ANY(dz(:) <= 0._wp)) CALL finish(routine, 'Dz must be larger than zero')
      new%dz     (:) = dz(:)
    END IF

    IF (PRESENT(lbounds) .AND. PRESENT(ubounds)) THEN
      IF (ANY(ubounds(:) < lbounds(:))) CALL finish(routine, 'Ubounds must be larger than Lbounds')
      IF (PRESENT(dz)) THEN
        IF (ANY(ABS((ubounds(:)-lbounds(:))-dz(:)) > EPSILON(1._wp))) CALL finish(routine, 'Dz != Ubounds-Lbounds')
      ELSE
        new%dz(:) = ubounds(:) - lbounds(:)
      END IF
      IF (.NOT. PRESENT(levels)) THEN
        new%levels(:) = 0.5_wp * (lbounds(:) + ubounds(:))
      END IF
    ELSE IF (PRESENT(lbounds)) THEN
      IF (PRESENT(dz)) THEN
        new%ubounds(:) = lbounds(:) + dz(:)
        IF (.NOT. PRESENT(levels)) THEN
          new%levels(:) = 0.5_wp * (lbounds(:) + new%ubounds(:))
        END IF
      ELSE IF (PRESENT(levels)) THEN
        new%ubounds(:) = levels(:) + (levels(:) - lbounds(:))
        new%dz(:) = new%ubounds(:) - lbounds(:)
      ELSE
        CALL finish(routine, 'Lbounds present ... specify either Levels or Dz!')
      END IF
    ELSE IF (PRESENT(ubounds)) THEN
      IF (PRESENT(dz)) THEN
        new%lbounds(:) = ubounds(:) - dz(:)
        IF (.NOT. PRESENT(levels)) THEN
          new%levels(:) = 0.5_wp * (new%lbounds(:) + ubounds(:))
        END IF
      ELSE IF (PRESENT(levels)) THEN
        new%lbounds(:) = levels(:) - (ubounds(:) - levels(:))
        new%dz(:) = ubounds(:) - new%lbounds(:)
      ELSE
        CALL finish(routine, 'Ubounds present ... specify either Levels or Dz!')
      END IF
    ELSE
      IF (PRESENT(levels) .AND. PRESENT(dz)) THEN
        CALL finish(routine, 'Levels and dz given ... please specify lbounds and/or ubounds')
      ELSE IF (PRESENT(dz)) THEN
        new%lbounds(1) = 0._wp
        DO i=2,length
          new%lbounds(i) = new%lbounds(i-1) + dz(i-1)
        END DO
        new%ubounds(:) = new%lbounds(:) + dz(:)
        new%levels(:) = 0.5_wp * (new%lbounds(:) + new%ubounds(:))
      ELSE IF (PRESENT(levels)) THEN
        IF (length > 1) THEN ! Multi-level axis
          new%lbounds(1) = 0._wp
          DO i=2,length
            new%lbounds(i) = levels(i-1) + (levels(i-1) - new%lbounds(i-1))
          END DO
          new%ubounds(:) = levels(:) + (levels(:) - new%lbounds(:))
          new%dz(:) = new%ubounds(:) - new%lbounds(:)
        ELSE ! Single level axis (level specified)
          new%lbounds(1) = levels(1)
          new%ubounds(1) = levels(1)
          new%dz(1)      = 0._wp
        END IF
      ELSE
        IF (length > 1) THEN ! Multi-level axis
          DO i=1,length
            new%levels(i) = REAL(i,wp)
            new%lbounds(i) = REAL(i,wp)
            new%ubounds(i) = REAL(i+1,wp)
            new%dz(i) = 1._wp
          END DO
        ELSE ! Single level axis (level not specified)
          new%levels(1) = 0._wp
          new%lbounds(1) = 0._wp
          new%ubounds(1) = 0._wp
          new%dz(1)      = 0._wp
        END IF
      END IF
    END IF
    !$ACC UPDATE DEVICE(new%levels) IF(present(levels))
    !$ACC UPDATE DEVICE(new%lbounds) IF(present(lbounds))
    !$ACC UPDATE DEVICE(new%ubounds) IF(present(ubounds))
    !$ACC UPDATE DEVICE(new%dz) IF(present(dz))
  END FUNCTION new_vgrid

  SUBROUTINE Register_vgrid(vgrid)

    USE mo_jsb_io_iface, ONLY: Create_zaxis

    TYPE(t_jsb_vgrid), POINTER, INTENT(inout) :: vgrid

    TYPE(t_jsb_vgrid_ptr), ALLOCATABLE :: vgrids_temp(:)
    INTEGER :: i, n, length

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
      n = SIZE(jsbach_vgrids)
      ALLOCATE(vgrids_temp(n+1))
      vgrids_temp(1:n) = jsbach_vgrids
      vgrids_temp(n+1) = t_jsb_vgrid_ptr(vgrid)
      CALL MOVE_ALLOC(vgrids_temp, jsbach_vgrids)
    ELSE
      ALLOCATE(jsbach_vgrids(1))
      jsbach_vgrids(1) = t_jsb_vgrid_ptr(vgrid)
    END IF
    vgrid%idx = SIZE(jsbach_vgrids)

    length = vgrid%n_levels
    CALL Create_zaxis(vgrid%ZaxisID, vgrid%echamZaxisIdx, &
      & vgrid%cdi_axis_type, length, vgrid%name, vgrid%levels(1:length), vgrid%units, vgrid%longname, &
      & vgrid%lbounds(1:length), vgrid%ubounds(1:length))

    WRITE(message_text, '(A,A,A,I10,I3)') 'Registered new vgrid: name, ID, length= ',TRIM(vgrid%name), ' ', vgrid%ZaxisID, length
    CALL message(TRIM(routine), message_text)

    !$ACC ENTER DATA CREATE(vgrid, vgrid%levels, vgrid%lbounds, vgrid%ubounds, vgrid%dz)
    !$ACC UPDATE DEVICE(vgrid%levels, vgrid%lbounds, vgrid%ubounds, vgrid%dz)
  END SUBROUTINE Register_vgrid

  FUNCTION Get_vgrid_by_name(name) RESULT(vgrid)

    CHARACTER(len=*),           INTENT(in)    :: name
    TYPE(t_jsb_vgrid), POINTER                :: vgrid

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Get_vgrid_by_name'

    vgrid => NULL()
    DO i=1,SIZE(jsbach_vgrids)
      IF (ASSOCIATED(jsbach_vgrids(i)%p)) THEN
        IF (jsbach_vgrids(i)%p%idx > 0 .AND. TRIM(jsbach_vgrids(i)%p%name) == TRIM(name)) THEN
          vgrid => jsbach_vgrids(i)%p
          EXIT
        END IF
      END IF
    END DO

    IF (.NOT. ASSOCIATED(vgrid)) THEN
      CALL finish(TRIM(routine), 'Vgrid '//TRIM(name)//' not found.')
    END IF

  END FUNCTION Get_vgrid_by_name

  FUNCTION Get_vgrid_by_idx_or_id(idx, id) RESULT(vgrid)

    INTEGER, INTENT(in), OPTIONAL :: idx, id
    TYPE(t_jsb_vgrid), POINTER :: vgrid

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_grid:Get_vgrid_by_idx_or_id'

    vgrid => NULL()

    IF (PRESENT(idx)) THEN
      IF (PRESENT(id)) CALL finish('routine', 'Only one of idx or id arguments allowed')
      IF (idx > SIZE(jsbach_vgrids)) CALL finish(routine, 'idx too large')
      IF (ASSOCIATED(jsbach_vgrids(idx)%p)) vgrid => jsbach_vgrids(idx)%p
    ELSE
      IF (PRESENT(idx)) CALL finish('routine', 'Only one of idx or id arguments allowed')
      IF (id >= INVALID) CALL finish(routine, 'vertical grid ID unknown')
      DO i=1,SIZE(jsbach_vgrids)
        IF (ASSOCIATED(jsbach_vgrids(i)%p)) THEN
          IF (jsbach_vgrids(i)%p%idx > 0 .AND. jsbach_vgrids(i)%p%id == id) THEN
            vgrid => jsbach_vgrids(i)%p
            EXIT
          END IF
        END IF
      END DO
    END IF
  
    IF (.NOT. ASSOCIATED(vgrid)) THEN
      CALL finish(TRIM(routine), 'Vgrid not found.')
    END IF

  END FUNCTION Get_vgrid_by_idx_or_id

  INTEGER FUNCTION Get_no_of_vgrids()

    Get_no_of_vgrids = SIZE(jsbach_vgrids)

  END FUNCTION Get_no_of_vgrids

  SUBROUTINE delete_vgrid(self)

    ! USE mo_jsb_io_iface, ONLY: Destroy_zaxis

#ifdef HAVE_F2003
    TYPE(t_jsb_vgrid), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_vgrid), POINTER :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':delete_vgrid'

    IF (self%id == 0) RETURN

    CALL message(routine, 'starting destruction of JSBACH vgrid instance.')

    !$ACC EXIT DATA DELETE(self%levels) IF(allocated(self%levels))
    !$ACC EXIT DATA DELETE(self%lbounds) IF(allocated(self%lbounds))
    !$ACC EXIT DATA DELETE(self%ubounds) IF(allocated(self%ubounds))
    !$ACC EXIT DATA DELETE(self%dz) IF(allocated(self%dz))
    IF (ALLOCATED(self%levels))  DEALLOCATE(self%levels)
    IF (ALLOCATED(self%lbounds)) DEALLOCATE(self%lbounds)
    IF (ALLOCATED(self%ubounds)) DEALLOCATE(self%ubounds)
    IF (ALLOCATED(self%dz))      DEALLOCATE(self%dz)

    DEALLOCATE(self)

    CALL message(routine, 'destruction of JSBACH vgrid completed.')

  END SUBROUTINE delete_vgrid

! Note: This subroutine is only called from ICON
  SUBROUTINE setup_zaxes_jsbach(vertical_axis_list)

    TYPE(t_verticalAxisList), INTENT(inout) :: vertical_axis_list

    TYPE(t_jsb_vgrid), POINTER :: vgrid
    INTEGER :: i

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::setup_zaxes_jsbach"

    IF (.NOT. ALLOCATED(jsbach_vgrids)) RETURN

    DO i=1,SIZE(jsbach_vgrids)
      vgrid => jsbach_vgrids(i)%p
      CALL Setup_jsb_vertical_axis(vertical_axis_list, &
        &                          vgrid%ZaxisID,      &
        &                          vgrid%name,         &
        &                          vgrid%n_levels,     &
        &                          vgrid%longname,     &
        &                          vgrid%units,        &
        &                          vgrid%levels,       &
        &                          vgrid%lbounds,      &
        &                          vgrid%ubounds       &
        &                         )
    END DO

  END SUBROUTINE setup_zaxes_jsbach

  SUBROUTINE set_vertical_grids_jsbach(vgrid_defs, count)

    TYPE(t_Vgrid), INTENT(inout) :: vgrid_defs(:)
    INTEGER,       INTENT(inout) :: count

    TYPE(t_jsb_vgrid), POINTER :: jsb_vgrid
    INTEGER :: i

    IF (.NOT. ALLOCATED(jsbach_vgrids)) RETURN

    DO i=1,SIZE(jsbach_vgrids)
      jsb_vgrid => jsbach_vgrids(i)%p
      CALL Set_jsb_restart_vgrid(vgrid_defs, count,  &
        &                        jsb_vgrid%ZaxisID,  &
        &                        jsb_vgrid%levels    &
        &                       )
    END DO

  END SUBROUTINE set_vertical_grids_jsbach

#endif
END MODULE mo_jsb_vertical_axes