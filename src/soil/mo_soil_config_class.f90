!>
!! @brief Contains structures and methods for soil config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_soil_config, new_soil_config

  TYPE t_soil_config
     LOGICAL :: active
     CHARACTER(len=filename_max) :: ic_filename
     CHARACTER(len=filename_max) :: bc_filename
     INTEGER :: nsoil_energy
     INTEGER :: nsoil_water
     REAL(wp), ALLOCATABLE :: soil_depth_energy(:)
     REAL(wp), ALLOCATABLE :: soil_depth_water(:)
     REAL(wp), ALLOCATABLE :: dz_energy(:)
     REAL(wp), ALLOCATABLE :: dz_water(:)
     REAL(wp) :: wsr_capacity               !< Water holding capacity of ground skin reservoir [m water equivalent]
     REAL(wp) :: wsn_capacity               !< Water holding capacity of ground snow reservoir (max snow depth) [m water equivalent]
  END type t_soil_config

  INTEGER, PARAMETER :: max_soil_layers = 20

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_config_class'

CONTAINS

  FUNCTION new_soil_config(nml_handler) RESULT(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    INTEGER,          INTENT(in) :: nml_handler
    TYPE(t_soil_config), POINTER :: config

    TYPE(t_soil_config), POINTER :: new

    LOGICAL  :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    INTEGER  :: nsoil_energy, nsoil_water
    REAL(wp) :: dz_energy(max_soil_layers)
    REAL(wp) :: dz_water(max_soil_layers)
    REAL(wp) :: wsr_capacity, wsn_capacity

    NAMELIST /jsb_soil_nml/                     &
      active,                                   &
      ic_filename, bc_filename,                 &
      nsoil_energy, nsoil_water,                &
      dz_energy, dz_water,                      &
      wsr_capacity, wsn_capacity

    INTEGER :: nml_unit, istat, i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_soil_config'

    CALL message(TRIM(routine), 'Starting soil configuration')

    ALLOCATE(new)

    ! Set defaults
    active         = .FALSE.
    ic_filename    = 'ic_land_soil.nc'
    bc_filename    = 'bc_land_soil.nc'
    nsoil_energy   = 5
    nsoil_water    = 5
    dz_energy(1:5) = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/) ; dz_energy(6:max_soil_layers) = 0._wp
    dz_water (1:5) = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/) ; dz_water (6:max_soil_layers) = 0._wp
    wsr_capacity   = 2.E-4_wp
    wsn_capacity   = 2.E-4_wp

    nml_unit = position_nml('jsb_soil_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_soil_nml)

    IF (nsoil_energy > max_soil_layers .OR. nsoil_water > max_soil_layers) THEN
      WRITE(message_text, '(A,I2)') 'Maximum number of soil layers is ', max_soil_layers
      CALL finish(TRIM(routine), message_text)
    END IF

    ALLOCATE(new%soil_depth_energy(nsoil_energy))
    ALLOCATE(new%soil_depth_water (nsoil_water ))
    ALLOCATE(new%dz_energy(nsoil_energy))
    ALLOCATE(new%dz_water (nsoil_water ))

    new%active                    = active
    new%ic_filename               = ic_filename
    new%bc_filename               = bc_filename
    new%nsoil_energy              = nsoil_energy
    new%nsoil_water               = nsoil_water
    new%soil_depth_energy(1)      = dz_energy(1)
    DO i=2,nsoil_energy
      new%soil_depth_energy(i) = new%soil_depth_energy(i-1) + dz_energy(i)
    END DO
    new%soil_depth_water(1)       = dz_water(1)
    DO i=2,nsoil_water
      new%soil_depth_water(i) = new%soil_depth_water(i-1) + dz_water(i)
    END DO
    new%dz_energy(1:nsoil_energy) = dz_energy(1:nsoil_energy)
    new%dz_water (1:nsoil_water)  = dz_water(1:nsoil_water)
    new%wsr_capacity              = wsr_capacity
    new%wsn_capacity              = wsn_capacity

    config => new

  END FUNCTION new_soil_config

#endif
END MODULE mo_soil_config_class
