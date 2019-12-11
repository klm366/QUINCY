!
! This file was automatically generated from JSBACH source file:
!    src/assimilation/mo_assimi_config_class.f90
! with md5 hash: c237bfc9cbdaa9147a82782d0f39aa87
!
!>
!! @brief Contains structures and methods for assimi config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (13.04.2016)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_assimi_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message !, message_text, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_assimi_config

  TYPE, EXTENDS(t_jsb_config) :: t_assimi_config
     INTEGER               :: ncanopy          ! number of canopy layers for the assimi process
                                               ! R: could be different for other processes e.g. for radiation
     REAL(wp), ALLOCATABLE :: canopy_bound_lai(:)
   CONTAINS
     PROCEDURE :: Init => Init_assimi_config
  END type t_assimi_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize assimi process
  !!
  !! @param[inout]     config     Configuration type of process (t_assimi_config)
  !!
  SUBROUTINE Init_assimi_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    CLASS(t_assimi_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    INTEGER                     :: ncanopy, i
    TYPE(t_jsb_vgrid), POINTER  :: canopy_layer

    NAMELIST /jsb_assimi_nml/  &
     & active,                 &
     & ic_filename,            &
     & bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_assimi_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting assimi configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_assimi.nc'
    ic_filename         = 'ic_land_assimi.nc'
    ncanopy             = 3

    !X e.g: use_alb_veg_simple        = .FALSE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_assimi_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_assimi_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%ncanopy             = ncanopy

    ! Create vertical axis
    canopy_layer  => new_vgrid('canopy_layer', ZAXIS_GENERIC, config%ncanopy, &
      levels=(/ (REAL(i,kind=wp),i=1,config%ncanopy) /), units='')
    CALL register_vgrid(canopy_layer)

    ! Prepare canopy boundaries of lai_cl
    ALLOCATE(config%canopy_bound_lai(0:ncanopy))
    DO i=0,ncanopy
       config%canopy_bound_lai(i) = REAL(i,wp) / REAL(ncanopy,wp)
       !write(*,*) "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR: canopy_bound_lai(i)" , new%canopy_bound_lai(i)
    END DO

  END SUBROUTINE Init_assimi_config

#endif
END MODULE mo_assimi_config_class
