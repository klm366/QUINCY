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
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_soil_config

  TYPE, EXTENDS(t_jsb_config) :: t_soil_config
     INTEGER :: nsoil_energy
     INTEGER :: nsoil_water
     REAL(wp), ALLOCATABLE :: soil_depth_energy(:)
     REAL(wp), ALLOCATABLE :: soil_depth_water(:)
     REAL(wp), ALLOCATABLE :: dz_energy(:)  !< layer thickness (dz means "delta z")
     REAL(wp), ALLOCATABLE :: dz_water(:)   !< layer thickness (dz means "delta z")
     REAL(wp) :: wsr_capacity               !< Water holding capacity of ground skin reservoir [m water equivalent]
     REAL(wp) :: wsn_capacity               !< Water holding capacity of ground snow reservoir (max snow depth) [m water equivalent]
     REAL(wp) :: ws_wilt_fract              !< Fractional soil moisture content at the wilting point (fraction of maximum moisture)
   CONTAINS
     PROCEDURE :: Init => Init_soil_config
  END type t_soil_config

  INTEGER, PARAMETER :: max_soil_layers = 20

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_config_class'

CONTAINS

  SUBROUTINE Init_soil_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND

    CLASS(t_soil_config), INTENT(inout) :: config

    LOGICAL  :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    INTEGER  :: nsoil_energy, nsoil_water
    REAL(wp) :: dz_energy(max_soil_layers)
    REAL(wp) :: dz_water(max_soil_layers)
    REAL(wp) :: wsr_capacity, wsn_capacity, ws_wilt_fract

    NAMELIST /jsb_soil_nml/                     &
      active,                                   &
      ic_filename, bc_filename,                 &
      nsoil_energy, nsoil_water,                &
      dz_energy, dz_water,                      &
      wsr_capacity, wsn_capacity, ws_wilt_fract

    INTEGER :: nml_handler, nml_unit, istat, i
    CHARACTER(len=filename_max) :: namelist_filename

    TYPE(t_jsb_vgrid), POINTER :: soil_e, soil_w

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_soil_config'

    CALL message(TRIM(routine), 'Starting soil configuration')

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
    ws_wilt_fract  = 0.35_wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_soil_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_soil_nml)

    CALL close_nml(nml_handler)

    IF (nsoil_energy > max_soil_layers .OR. nsoil_water > max_soil_layers) THEN
      WRITE(message_text, '(A,I2)') 'Maximum number of soil layers is ', max_soil_layers
      CALL finish(TRIM(routine), message_text)
    END IF

    ALLOCATE(config%soil_depth_energy(nsoil_energy))
    ALLOCATE(config%soil_depth_water (nsoil_water ))
    ALLOCATE(config%dz_energy(nsoil_energy))
    ALLOCATE(config%dz_water (nsoil_water ))

    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%nsoil_energy              = nsoil_energy
    config%nsoil_water               = nsoil_water
    config%soil_depth_energy(1)      = dz_energy(1)
    DO i=2,nsoil_energy
      config%soil_depth_energy(i) = config%soil_depth_energy(i-1) + dz_energy(i)
    END DO
    config%soil_depth_water(1)       = dz_water(1)
    DO i=2,nsoil_water ! R: Note: if nsoil_water is 1 this loop will be skiped automatically
      config%soil_depth_water(i) = config%soil_depth_water(i-1) + dz_water(i)
    END DO
    config%dz_energy(1:nsoil_energy) = dz_energy(1:nsoil_energy)
    config%dz_water (1:nsoil_water)  = dz_water(1:nsoil_water)
    config%wsr_capacity              = wsr_capacity
    config%wsn_capacity              = wsn_capacity
    config%ws_wilt_fract             = ws_wilt_fract
    IF (.NOT. active) RETURN

    ! Create vertical axes

    soil_e  => new_vgrid('soil_depth_energy', ZAXIS_DEPTH_BELOW_LAND, config%nsoil_energy, &
      levels=config%soil_depth_energy(1:config%nsoil_energy), units='m')
    CALL register_vgrid(soil_e)

    soil_w  => new_vgrid('soil_depth_water', ZAXIS_DEPTH_BELOW_LAND, config%nsoil_water, &
      levels=config%soil_depth_water(1:config%nsoil_water), units='m')
    CALL register_vgrid(soil_w)

  END SUBROUTINE Init_soil_config

#endif
END MODULE mo_soil_config_class
