! from the JSBACH4 git develop branch rev. 972a425 - Steffen 01.09.2017

!>
!! @brief Contains JSBACH grid structure and methods
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_grid_class
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message_text, message
  USE mo_jsb_domain_iface,   ONLY: t_patch, pe_decomposed

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_grid, t_jsb_vgrid

  TYPE t_jsb_grid
    INTEGER                      :: id = 0
    CHARACTER(LEN=10)            :: name = ''
    CHARACTER(LEN=filename_max)  :: filename = ''
    INTEGER                      :: ntotal = 0          !< Number of cells on this PE
    INTEGER                      :: nland = 0
    REAL(wp), POINTER            :: lon(:,:) => NULL()
    REAL(wp), POINTER            :: lat(:,:) => NULL()
    REAL(wp), POINTER            :: area(:,:) => NULL()
    LOGICAL,  POINTER            :: lsm(:,:) => NULL()
    REAL(wp), POINTER            :: lsf(:,:) => NULL()
    INTEGER                      :: nblks = 1
    INTEGER                      :: nproma = 0
    INTEGER                      :: npromz = 0
    INTEGER                      :: ntotal_g              !< Total number of cells in grid
    INTEGER                      :: dims_g(2)             !< Dimensions of global grid
    INTEGER                      :: nlat_g = 0            !< Number of latitudes in global grid (effective for ICON)
    LOGICAL                      :: configured = .FALSE.
    TYPE(t_patch), POINTER       :: patch
    TYPE(pe_decomposed), POINTER :: decomposition
  CONTAINS
    PROCEDURE :: Print => Print_grid
  END TYPE t_jsb_grid

  TYPE t_jsb_vgrid
    INTEGER ::           id           = 0
    CHARACTER(len=20) :: name         = ''
    CHARACTER(len=50) :: longname     = ''
    CHARACTER(len=10) :: units        = ''
    INTEGER :: n_levels               = 0
    REAL(wp), ALLOCATABLE :: levels(:)
    REAL(wp), ALLOCATABLE :: lbounds(:)
    REAL(wp), ALLOCATABLE :: ubounds(:)
    REAL(wp), ALLOCATABLE :: dz(:)
    INTEGER :: cdiZaxisID             = 0
    INTEGER :: echamZaxisID           = 0
  END TYPE t_jsb_vgrid

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_grid_class'

CONTAINS

  !>
  !! Print method for grid instance.
  !!
  !! Subroutine to print description of grid instance.
  !!
  SUBROUTINE Print_grid(this)

    CLASS(t_jsb_grid), INTENT(in) :: this

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Print_grid'

    CALL message(TRIM(routine), 'print JSBACH grid description')

    WRITE(message_text,'(A,I2,A,A)') 'Grid id: ', this%id, ', grid name: ', TRIM(this%name)
    CALL message('  ', message_text)

    IF (.NOT. this%configured) THEN
      CALL message(TRIM(routine), 'Grid not configured yet!')
      RETURN
    END IF

    IF (this%filename /= '') THEN
      WRITE(message_text,'(A)') 'Grid filename: '//TRIM(this%filename)
      CALL message('  ', message_text)
    END IF

    IF (this%patch%id >= 0) THEN
      WRITE(message_text,'(A,I2)') 'Patch ID: ', this%patch%id
      CALL message('  ', message_text)
    END IF

  END SUBROUTINE Print_grid

#endif
END MODULE mo_jsb_grid_class
