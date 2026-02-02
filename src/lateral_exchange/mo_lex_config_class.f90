!> Contains structures and methods for lateral exchange process configuration
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Configuration of the lateral exchange process
!>
!> The module defines the data structure [[t_lex_config]] that contains the lateral exchange
!> (LEX) configuration, i.e. all lateral exchange variables that are configurable via namelist.
!> In subroutine [[Init_lex_config]] we read the lateral exchange namelist 'jsb_lex_nml'
!> and set the lateral exchange configuration data structure accordingly. We also perform
!> some consistency checks.
!>
MODULE mo_lex_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_lex_config
  !> The data structure contains all configurable parameters of the LEX process. It extends
  !> [[mo_jsb_config_class:t_jsb_config]], a data structure containing basic parameters defined
  !> for all ICON-Land processes.
  !>
  !> _Note for the ford documentation_: The table lists the general parameters of
  !> [[mo_jsb_config_class:t_jsb_config]] as well as the parameters specific to the lateral exchange
  !> actually defined in [[t_lex_config]].

  ! TODO: Currently a lot of switches are defined to conveniently turn on/off different lateral exchange
  !       processes for testing. In future versions of the model the amount of switches should be reduced.
  TYPE, EXTENDS(t_jsb_config) :: t_lex_config
    LOGICAL                   :: &
      & l_runoff_trns,           & !< en/disable lateral surface runoff to connected tiles
      & l_drainage_trns,         & !< en/disable lateral drainage to connected tiles
      & l_infil_trns,            & !< en/disable lateral infiltration of pond water in neighbouring soil tile
      & l_snow_trns,             & !< en/disable lateral transport of snow between tiles
      & l_dffsn_wtr,             & !< en/disable diffusive water exchange between tiles
      & l_dffsn_nrg,             & !< en/disable diffusive heat exchange between tiles
      & l_cnvctn_nrg,            & !< en/disable advective (via movement of water) heat exchange between tiles
      & l_wellmixed,             & !< T: tile fractions are assumed to be well-mixed in the grid cell
      & l_prtlcon,               & !< T: tiles can be only partially connected to downstream tiles
      & l_macrtrns,              & !< en/disable macro-scale lateral transport between tiles
      & l_micrtrns,              & !< en/disable micro-scale lateral transport between tiles
      & l_nested,                & !< T: nested setup in which a tile can be surrounded by another tile, which
                                         again can be surrounded by other tiles
      & l_dpr4ponds,             & !< T: depression characteristics are used for pond characteristics
                                   !< if no tile-specific pond fraction / depth is provided
      & l_noslope_dpr              !< T: if no tile-specific macro-slope is provided, the slope on
                                   !< depression tiles (h_surf < 0) is set to zero
    INTEGER                   :: &
      & ncontis                    !< Number of connected tiles used in setup

   CONTAINS
     PROCEDURE :: Init => Init_lex_config
  END type t_lex_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lex_config_class'

CONTAINS
  ! ============================================================================================== !
  !>
  !>#### Configuration of the lateral exchange process
  !>
  !> In this subroutine we configure the lateral exchange process based on the settings in the lateral
  !> exchange namelist 'jsb_lex_nml'. Init_lex_config is a type-bound procedure of [[t_lex_config]].
  !> It is called during the initialization phase from [[mo_jsb_model_class:Configure_model_processes]]
  !> ('CALL this%processes(iproc)%p%Configure()').
  !>
  !> _Note_: For a more extensive documentation of the namelist parameters please check the
  !> configuration data type [[t_lex_config]].
  !>
  SUBROUTINE Init_lex_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    CLASS(t_lex_config), INTENT(inout) :: config

    !--------------------------------------------------------------------------------------------------------
    ! Namelist parameters
    !
    LOGICAL                     :: &
      & active                      !< T: LEX process is active
    INTEGER                     :: &
      & ncontis                      !< Number of connected tiles
    CHARACTER(len=filename_max) :: &
      & ic_filename,               & !< Name of input file with initial conditions
      & bc_filename,               & !< Name of input file with boundary conditions
    LOGICAL                     :: &
      & l_runoff_trns,           & !< en/disable lateral surface runoff to connected tiles
      & l_drainage_trns,         & !< en/disable lateral drainage to connected tiles
      & l_infil_trns,            & !< en/disable lateral infiltration of pond water in neighbouring soil tile
      & l_snow_trns,             & !< en/disable lateral transport of snow between tiles
      & l_dffsn_wtr,             & !< en/disable diffusive water exchange between tiles
      & l_dffsn_nrg,             & !< en/disable diffusive heat exchange between tiles
      & l_cnvctn_nrg,            & !< en/disable advective (via movement of water) heat exchange between tiles
      & l_wellmixed,             & !< T: tile fractions are assumed to be well-mixed in the grid cell
      & l_prtlcon,               & !< T: tiles can be only partially connected to downstream tiles
      & l_macrtrns,              & !< en/disable macro-scale lateral transport between tiles
      & l_micrtrns,              & !< en/disable micro-scale lateral transport between tiles
      & l_nested,                & !< T: nested setup in which a tile can be surrounded by another tile,
                                   !< which again my be surrounded by other tiles
      & l_dpr4ponds,             & !< T: depression characteristics are used for pond characteristics
                                   !< if no tile-specific pond fraction / depth is provided
      & l_noslope_dpr              !< T: if no tile-specific macro-slope is provided, the slope on
                                   !< depression tiles (h_surf < 0) is set to zero

    NAMELIST /jsb_lex_nml/      &
         active,                &
         ncontis,               &
         l_runoff_trns,         &
         l_drainage_trns,       &
         l_infil_trns,          &
         l_snow_trns,           &
         l_dffsn_wtr,           &
         l_dffsn_nrg,           &
         l_cnvctn_nrg,          &
         l_wellmixed,           &
         l_prtlcon,             &
         l_macrtrns,            &
         l_micrtrns,            &
         l_nested,              &
         l_dpr4ponds,           &
         l_noslope_dpr

    INTEGER ::       &
      & nml_handler, &  !< Handle to deal with namelist file
      & nml_unit,    &  !< IO unit of the namelist file
      & istat,       &  !< Status flag
      & i               !< Loop variable

    TYPE(t_jsb_vgrid), POINTER  :: tile_layer !< Grid axis used for connected tiles

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lex_config'

    CALL message(TRIM(routine), 'Starting lex configuration')

    !>
    !> Default values are defined for all namelist parameters.
    !
    active              = .FALSE.  !< LEX process is inactive by default
    ncontis             = 0
    l_runoff_trns       = .TRUE.
    l_drainage_trns     = .TRUE.
    l_infil_trns        = .TRUE.
    l_snow_trns         = .TRUE.
    l_dffsn_wtr         = .TRUE.
    l_dffsn_nrg         = .TRUE.
    l_cnvctn_nrg        = .TRUE.
    l_wellmixed         = .TRUE.
    l_prtlcon           = .TRUE.
    l_macrtrns          = .TRUE.
    l_micrtrns          = .TRUE.
    l_nested            = .TRUE.
    l_dpr4ponds         = .TRUE.
    l_noslope_dpr       = .TRUE.

    ! Open the file containing the lateral exchange namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    ! Find and read the lateral exchange namelist
    ! The parameters defined in the namelist overwrite the default values defined above.
    nml_unit = position_nml('jsb_lex_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_lex_nml)

    CALL close_nml(nml_handler)

    !> After reading the namelist, all configurable parameters for the lateral exchange - default
    !> parameters and values set in the namelist - are stored in the configuration data structure
    !> 'config' to make them easily accessible in all ICON-Land process modules ("dsl4jsb_Use_config(LEX_)").
    !>
    !> The specific choice of parameters is printed to the experiment's log file to document the
    !> setup. The simulation is stopped if the number of connected tiles is not provided via namelist.

    config%active              = active
    config%l_runoff_trns       = l_runoff_trns
    config%l_drainage_trns     = l_drainage_trns
    config%l_infil_trns        = l_infil_trns
    config%l_snow_trns         = l_snow_trns
    config%l_dffsn_wtr         = l_dffsn_wtr
    config%l_dffsn_nrg         = l_dffsn_nrg
    config%l_cnvctn_nrg        = l_cnvctn_nrg
    config%l_wellmixed         = l_wellmixed
    config%l_prtlcon           = l_prtlcon
    config%l_macrtrns          = l_macrtrns
    config%l_micrtrns          = l_micrtrns
    config%l_nested            = l_nested
    config%l_dpr4ponds         = l_dpr4ponds
    config%l_noslope_dpr       = l_noslope_dpr

    IF (.NOT. active) RETURN

    IF (ncontis == 0) THEN
      CALL finish(TRIM(routine), 'Number of connected tiles needs to be defined via namelist.')
    ELSE
      config%ncontis = ncontis
    END IF

    IF (l_runoff_trns) THEN
      CALL message(TRIM(routine), 'Lateral surface runoff is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral surface runoff is disabled')
    END IF

    IF (l_drainage_trns) THEN
      CALL message(TRIM(routine), 'Lateral drainage is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral drainage is disabled')
    END IF

    IF (l_infil_trns) THEN
      CALL message(TRIM(routine), 'Lateral infiltration is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral infiltration is disabled')
    END IF

    IF (l_snow_trns) THEN
      CALL message(TRIM(routine), 'Lateral snow redistribution between tiles is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral snow redistribution between tiles is disabled')
    END IF

    IF (l_dffsn_wtr) THEN
      CALL message(TRIM(routine), 'Lateral diffusion (water) is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral diffusion (water) is disabled')
    END IF

    IF (l_dffsn_nrg) THEN
      CALL message(TRIM(routine), 'Lateral diffusion (energy) is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral diffusion (energy) is disabled')
    END IF

    IF (l_cnvctn_nrg) THEN
      CALL message(TRIM(routine), 'Lateral convection (energy) is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral convection (energy) is disabled')
    END IF

    IF (l_wellmixed) THEN
      CALL message(TRIM(routine), 'Tile fractions are assumed to be well-mixed within the grid cell')
    ELSE
      CALL message(TRIM(routine), 'Microcon_ctc based exclusivly on values of matrix elements')
    END IF

    IF (l_prtlcon) THEN
      CALL message(TRIM(routine), 'Allowing for tiles to be only partially connected to downstream tiles')
    ELSE
      CALL message(TRIM(routine), 'Tiles are either fully or not at all connected to downstream tiles')
    END IF

    IF (l_nested) THEN
      CALL message(TRIM(routine), 'Assuming a nested setup in which an encompassed tile may sourround another tile')
    ELSE
      CALL message(TRIM(routine), 'Nested constellations are being ignored and '// &
        & 'geometry is calculated exclusivly based on the characteristics of two tiles)')
    END IF

    IF (l_dpr4ponds) THEN
      CALL message(TRIM(routine), 'Depressions constitute ponding area (if no tile specific values are provided)')
    ELSE
      CALL message(TRIM(routine), 'Uniform ponding area is assumed for tiles (if no tile specific values are provided)')
    END IF

    IF (l_noslope_dpr) THEN
      CALL message(TRIM(routine), 'Zero slope is assumed for depressions (if no tile specific values are provided)')
    ELSE
      CALL message(TRIM(routine), 'Slope is set to steepness for depressions (if no tile specific values are provided)')
    END IF

    !< Create additional axis for connected tiles
    ! This axis is not a vertical axis in the usual sense. It is required to organize lateral interactions between tiles
    tile_layer  => new_vgrid('connecttiles', ZAXIS_GENERIC, ncontis, &
      & levels=(/ (REAL(i,kind=wp),i=1,ncontis) /), units='')
    CALL register_vgrid(tile_layer)

  END SUBROUTINE Init_lex_config

#endif
END MODULE mo_lex_config_class
