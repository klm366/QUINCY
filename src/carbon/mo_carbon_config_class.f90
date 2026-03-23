!> Contains structures and methods for carbon config
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_carbon_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message, message_text
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_carbon_config

  TYPE, EXTENDS(t_jsb_config) :: t_carbon_config
    LOGICAL  :: with_nitrogen
    LOGICAL  :: read_cpools ! Read Cpools from external file when experiment is newly started (cpools are
                            ! initialized with values from file) or restarted (cpools from restart file are
                            ! overwritten by values from file).
                            ! ATTENTION: The run script has to take care that read_cpools is set to .FALSE.
                            !            in the experiment jobs following the first one! Or re-set it manually!
    LOGICAL  :: diag_humus_fluxes ! Activate diagnostic carbon fluxes to and from the yasso humus pools
                                  ! required for analytical humus pool equilibration (compare jsb3.2 svn rev 9934)
    LOGICAL  :: on_soillev  ! T: Carbon calculated on vertical grid (soil layers)
    LOGICAL  :: clear_frozen
    REAL(wp) :: resp_rate_modifier
    REAL(wp) :: fire_fract_wood_2_atmos

   CONTAINS
     PROCEDURE :: Init => Init_carbon_config
  END type t_carbon_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize carbon process
  !!
  !! @param[inout]     config     Configuration type of process (t_carbon_config)
  !!
  SUBROUTINE Init_carbon_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_carbon_constants,   ONLY: FireFracWood2Atmos

    CLASS(t_carbon_config), INTENT(inout) :: config

    LOGICAL                     :: active, lrestart_cont
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: with_nitrogen
    LOGICAL                     :: read_cpools, diag_humus_fluxes
    LOGICAL                     :: on_soillev
    LOGICAL                     :: clear_frozen
    REAL(wp)                    :: resp_rate_modifier

    REAL(wp)                    :: fire_fract_wood_2_atmos
    ! LOGICAL           :: lemissions

    NAMELIST /jsb_carbon_nml/      &
         active,                   &
         lrestart_cont,            &
         ic_filename,              &
         bc_filename,              &
         with_nitrogen,            &
         read_cpools,              &
         diag_humus_fluxes,        &
         on_soillev,               &
         clear_frozen,             &
         resp_rate_modifier

    INTEGER :: nml_handler, nml_unit, istat
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)
    REAL(wp), ALLOCATABLE :: depths(:), mids(:)
    INTEGER :: nsoil

    TYPE(t_jsb_vgrid), POINTER :: soil_c

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_carbon_config'

    CALL message(TRIM(routine), 'Starting carbon configuration')

    ! Set defaults
    active                 = .FALSE.
    lrestart_cont          = .FALSE.       ! TRUE: Continue although CARBON variables are missing in restart file
    bc_filename            = 'bc_land_carbon.nc'
    ic_filename            = 'ic_land_carbon.nc'
    with_nitrogen          = .FALSE.
    read_cpools            = .FALSE.
    diag_humus_fluxes      = .FALSE.
    on_soillev             = .FALSE.
    clear_frozen           = .FALSE.
    resp_rate_modifier      = 1._wp
    fire_fract_wood_2_atmos = FireFracWood2Atmos

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_carbon_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_carbon_nml)

    CALL close_nml(nml_handler)

    ! Write resulting values in config memory
    config%active                 = active
    config%lrestart_cont          = lrestart_cont
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%with_nitrogen          = with_nitrogen
    config%read_cpools            = read_cpools
    config%diag_humus_fluxes      = diag_humus_fluxes
    config%on_soillev             = on_soillev
    config%clear_frozen           = clear_frozen
    config%resp_rate_modifier     = resp_rate_modifier
    config%fire_fract_wood_2_atmos = fire_fract_wood_2_atmos

    IF (on_soillev) THEN
      input_file = jsb_netcdf_open_input(ic_filename)

      ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
      nsoil = SIZE(ptr_1D)

      ALLOCATE(depths(nsoil+1))
      ALLOCATE(mids(nsoil))
      depths(1) = 0._wp
      depths(2:nsoil+1) = ptr_1D(1:nsoil)
      mids(1:nsoil) = (depths(1:nsoil) + depths(2:nsoil+1)) / 2._wp
      DEALLOCATE(ptr_1D)

      CALL input_file%Close()

      soil_c => new_vgrid('soil_depth_carbon', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
        & levels  = mids                 (1:nsoil  ),                          &
        & lbounds = depths               (1:nsoil  ),                          &
        & ubounds = depths               (2:nsoil+1),                          &
        & units='m')

      CALL register_vgrid(soil_c)
      WRITE(message_text, *) 'Soil levels in carbon (upper) [m]: ', soil_c%lbounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Soil levels in carbon (mid)   [m]: ', soil_c%levels
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Soil levels in carbon (lower) [m]: ', soil_c%ubounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Soil level depths in carbon   [m]: ', soil_c%dz
      CALL message(TRIM(routine), message_text)
      DEALLOCATE(depths, mids)
    END IF

  END SUBROUTINE Init_carbon_config

#endif
END MODULE mo_carbon_config_class
