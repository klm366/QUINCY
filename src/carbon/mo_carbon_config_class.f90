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
!>#### Configuration of the carbon process
!>
!> The module defines the data structure [[t_carbon_config]] that contains the carbon process
!> configuration, i.e. all carbon variables that are configurable via namelist.
!> In subroutine [[Init_carbon_config]] we read the carbon namelist 'jsb_carbon_nml'
!> and set the carbon configuration data structure accordingly. We also perform
!> some consistency checks.
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
  !> The data structure contains all configurable parameters of the CARBON process. It extends
  !> [[mo_jsb_config_class:t_jsb_config]], a data structure containing basic parameters defined
  !> for all ICON-Land processes.
  !>
  !> _Note for the ford documentation_: The table lists the general parameters of
  !> [[mo_jsb_config_class:t_jsb_config]] as well as the parameters specific to the carbon process
  !> actually defined in [[t_carbon_config]].

  TYPE, EXTENDS(t_jsb_config) :: t_carbon_config
    LOGICAL  :: with_nitrogen
    LOGICAL  :: read_cpools !< Read Cpools from external file when experiment is newly started (cpools are
                            !< initialized with values from file) or restarted (cpools from restart file are
                            !< overwritten by values from file).
                            !< ATTENTION: The run script has to take care that read_cpools is set to .FALSE.
                            !<            in the experiment jobs following the first one! Or re-set it manually!
    LOGICAL  :: diag_humus_fluxes  !< Activate diagnostic carbon fluxes to and from the yasso humus pools
                                   !< required for analytical humus pool equilibration (compare jsb3.2 svn rev 9934)
    LOGICAL  :: vertical_yasso     !< T: Below ground carbon calculated on vertical grid (soil layers)
    LOGICAL  :: remove_permafrost_carbon !< T: Remove carbon from permafrost affected soil layers at simulation (re)start
    REAL(wp) :: fire_fract_wood_2_atmos !< Fraction of wood carbon emitted to the atm by fire.

   CONTAINS
     PROCEDURE :: Init => Init_carbon_config
  END type t_carbon_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_config_class'

CONTAINS

  ! ============================================================================================== !
  !>
  !>#### Configuration of the carbon process
  !>
  !> In this subroutine we configure the carbon process based on the settings in the carbon
  !> namelist 'jsb_carbon_nml'. Init_carbon_config is a type-bound procedure of [[t_carbon_config]].
  !> It is called during initialization phase from [[mo_jsb_model_class:Configure_model_processes]]
  !> ('CALL this%processes(iproc)%p%Configure()').
  !>
  !> _Note_: For a more extensive documentation of the namelist parameters please check the
  !> configuration data type [[t_carbon_config]].
  !>
  SUBROUTINE Init_carbon_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_carbon_constants,   ONLY: FireFracWood2Atmos

    CLASS(t_carbon_config), INTENT(inout) :: config

    LOGICAL ::                     &
      & active,                    & !< True: the CARBON process is activated
      & lrestart_cont                !< True: continue simulations even if CARBON restart data is missing
    CHARACTER(len=filename_max) :: &
      & ic_filename,               & !< Name of input file with initial conditions
      & bc_filename                  !< Name of input file with boundary conditions
    LOGICAL ::                     &
      & with_nitrogen,             &
      & read_cpools,               & !< Read Cpools from external file when experiment is newly started (cpools are
                                     !< initialized with values from file) or restarted (cpools from restart file are
                                     !< overwritten by values from file).
                                     !< ATTENTION: The run script has to take care that read_cpools is set to .FALSE.
                                     !<            in the experiment jobs following the first one! Or re-set it manually!
      & diag_humus_fluxes,         & !< Activate diagnostic carbon fluxes to and from the yasso humus pools
                                     !< required for analytical humus pool equilibration (compare jsb3.2 svn rev 9934)
      & vertical_yasso,            & !< T: Below ground carbon calculated on vertical grid (soil layers)
      & remove_permafrost_carbon     !< T: Remove carbon from permafrost affected soil layers at simulation (re)start
    REAL(wp) ::                    &
      & fire_fract_wood_2_atmos      !< Fraction of wood carbon emitted to the atm by fire.

    NAMELIST /jsb_carbon_nml/      &
         active,                   &
         lrestart_cont,            &
         ic_filename,              &
         bc_filename,              &
         with_nitrogen,            &
         read_cpools,              &
         diag_humus_fluxes,        &
         vertical_yasso,           &
         remove_permafrost_carbon

    INTEGER ::                     &
      & nml_handler,               &   !< Handle to deal with namelist file
      & nml_unit,                  &   !< IO unit of the namelist file
      & istat                          !< Status flag
    TYPE(t_input_file) :: input_file   !< Input file
    REAL(wp), POINTER :: ptr_1D(:)     !< Pointer to a one dimensional variable
    REAL(wp), ALLOCATABLE :: depths(:) !< Upper/lower bounds of soil layers
    REAL(wp), ALLOCATABLE :: mids(:)   !< Mid soil layer depth
    INTEGER :: nsoil                   !< Number of soil layers

    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_c !< Vertical grid used for soil layers of CARBON process

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_carbon_config'

    CALL message(TRIM(routine), 'Starting carbon configuration')

    !>
    !> Default values are defined for all namelist parameters.
    !
    active                   = .FALSE.
    lrestart_cont            = .FALSE.       ! TRUE: Continue although CARBON variables are missing in restart file
    bc_filename              = 'bc_land_carbon.nc'
    ic_filename              = 'ic_land_carbon.nc'
    with_nitrogen            = .FALSE.
    read_cpools              = .FALSE.
    diag_humus_fluxes        = .FALSE.
    vertical_yasso           = .FALSE.
    remove_permafrost_carbon = .FALSE.
    fire_fract_wood_2_atmos  = FireFracWood2Atmos

    ! Open the file containing the carbon namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    ! Find and read the carbon namelist
    ! The parameters defined in the namelist overwrite the default values defined above.
    nml_unit = position_nml('jsb_carbon_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_carbon_nml)

    CALL close_nml(nml_handler)

    !> After reading the namelist, all configurable parameters of the carbon - default parameters
    !> and values set in the namelist - are stored in the configuration data structure 'config' to
    !> make them easily accessible in all ICON-Land process modules ("dsl4jsb_Use_config(CARBON_)").

    config%active                   = active
    config%lrestart_cont            = lrestart_cont
    config%ic_filename              = ic_filename
    config%bc_filename              = bc_filename
    config%with_nitrogen            = with_nitrogen
    config%read_cpools              = read_cpools
    config%diag_humus_fluxes        = diag_humus_fluxes
    config%vertical_yasso           = vertical_yasso
    config%remove_permafrost_carbon = remove_permafrost_carbon
    config%fire_fract_wood_2_atmos  = fire_fract_wood_2_atmos

    IF (vertical_yasso) THEN
      input_file = jsb_netcdf_open_input(ic_filename)

      ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
      nsoil = SIZE(ptr_1D)                                     ! Number of soil layers

      ! calculate mid-points of each soil layer
      ALLOCATE(depths(nsoil+1))                                     ! Upper/lower bounds of soil layers
      ALLOCATE(mids(nsoil))
      depths(1) = 0._wp                                             ! lower bounds == top of the soil layer
      depths(2:nsoil+1) = ptr_1D(1:nsoil)                           ! upper bounds == bottom of the soil layer
      mids(1:nsoil) = (depths(1:nsoil) + depths(2:nsoil+1)) / 2._wp ! center depth of a layer
      DEALLOCATE(ptr_1D)

      CALL input_file%Close()

      vgrid_soil_c => new_vgrid('soil_depth_carbon', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
        & levels  = mids                 (1:nsoil  ),                          &
        & lbounds = depths               (1:nsoil  ),                          &
        & ubounds = depths               (2:nsoil+1),                          &
        & units='m')
      CALL register_vgrid(vgrid_soil_c)

      DEALLOCATE(depths, mids)

      CALL message(TRIM(routine), 'Init vgrid soil_carbon_energy with values from input file: '//TRIM(ic_filename))

      ! Info saved to logfile
      WRITE(message_text, *) 'Soil levels in carbon (upper) [m]: ', vgrid_soil_c%lbounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Soil levels in carbon (mid)   [m]: ', vgrid_soil_c%levels
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Soil levels in carbon (lower) [m]: ', vgrid_soil_c%ubounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Soil level depths in carbon   [m]: ', vgrid_soil_c%dz
      CALL message(TRIM(routine), message_text)
    END IF

  END SUBROUTINE Init_carbon_config

#endif
END MODULE mo_carbon_config_class
