!> Contains structures and methods for lex config
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
!>#### Contains namelist information of lex
!>
MODULE mo_lex_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_lex_config

  TYPE, EXTENDS(t_jsb_config) :: t_lex_config
     LOGICAL                  :: &
       & l_runoff_trns,          & !< enable lateral surface runoff between tiles
       & l_drainage_trns,        & !< enable lateral drainage between tiles
       & l_infil_trns,           & !< enable lateral infiltration between tiles
       & l_snow_trns,            & !< enable lateral transport of snow between tiles
       & l_dffsn_wtr,            & !< enable diffusive water exchange between tiles
       & l_dffsn_nrg,            & !< enable diffusive heat exchange between tiles
       & l_nested,               & !< Assuming a nested constellation in which an encompassed tile may sourround 1 tile itself
       & l_dpr4ponds,            & !< if no tile-specific pond fraction / depth is provided assume depressions to be ponds  
       & l_noslope_dpr             !< if no tile-specific macro-slope is provided, the slope on depression tiles (h_surf < 0) is set to zero 
    INTEGER                     :: ncontis !< number of connected tiles used in setup

   CONTAINS
     PROCEDURE :: Init => Init_lex_config
  END type t_lex_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lex_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize lex process
  !!
  SUBROUTINE Init_lex_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    CLASS(t_lex_config), INTENT(inout) :: config

    LOGICAL                     :: active
    INTEGER                     :: ncontis
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: &
       & l_runoff_trns,            & !< enable lateral surface runoff between tiles
       & l_drainage_trns,          & !< enable lateral drainage between tiles
       & l_infil_trns,             & !< enable lateral infiltration between tiles
       & l_snow_trns,              & !< enable lateral transport of snow between tiles
       & l_dffsn_wtr,              & !< enable diffusive water exchange between tiles
       & l_dffsn_nrg,              & !< enable diffusive heat exchange between tiles
       & l_nested,                 & !< Assuming a nested constellation in which an encompassed tile may sourround 1 tile itself 
       & l_dpr4ponds,              & !< if no tile -specific pond fraction / depth is provided assume depressions to be ponds 
       & l_noslope_dpr               !< if no tile -specific macro-slope is provided, the slope on depression tiles (h_surf < 0) is set to zero  

    TYPE(t_jsb_vgrid), POINTER  :: tile_layer             

    NAMELIST /jsb_lex_nml/      &
         active,                & !< process is active
         ncontis,               & !< number of connected tiles used in setup
         l_runoff_trns,         & !< enable lateral surface runoff between tiles
         l_drainage_trns,       & !< enable lateral drainage between tiles
         l_infil_trns,          & !< enable lateral infiltration between tiles
         l_snow_trns,           & !< enable lateral transport of snow between tiles
         l_dffsn_wtr,           & !< enable diffusive water exchange between tiles
         l_dffsn_nrg,           & !< enable diffusive heat exchange between tiles
         l_nested,              & !< Assuming a nested constellation in which an encompassed tile may sourround 1 tile itself 
         l_dpr4ponds,           & !< if no tile -specific pond fraction / depth is provided assume depressions to be ponds 
         l_noslope_dpr            !< if no tile -specific macro-slope is provided, the slope on depression tiles (h_surf < 0) is set to zero   

    INTEGER :: nml_handler, nml_unit, istat, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lex_config'

    CALL message(TRIM(routine), 'Starting lex configuration')


    ! Set defaults
    active              = .FALSE.
    ncontis             = 0
    l_runoff_trns       = .TRUE.
    l_drainage_trns     = .TRUE.
    l_infil_trns        = .TRUE.
    l_snow_trns         = .TRUE.
    l_dffsn_wtr         = .TRUE.
    l_dffsn_nrg         = .TRUE.
    l_nested            = .TRUE.
    l_dpr4ponds         = .TRUE.
    l_noslope_dpr       = .TRUE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_lex_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_lex_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%l_runoff_trns       = l_runoff_trns
    config%l_drainage_trns     = l_drainage_trns
    config%l_infil_trns        = l_infil_trns
    config%l_snow_trns         = l_snow_trns
    config%l_dffsn_wtr         = l_dffsn_wtr
    config%l_dffsn_nrg         = l_dffsn_nrg
    config%l_nested            = l_nested
    config%l_dpr4ponds         = l_dpr4ponds
    config%l_noslope_dpr       = l_noslope_dpr

    IF (.NOT. active) RETURN

    ! get number of connected tiles from model config
    IF (ncontis == 0) THEN
      CALL finish(TRIM(routine), 'Trying to configure lateral exchange between tiles, without connected tiles.')
    ELSE
      config%ncontis           = ncontis  
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

    IF (l_nested) THEN
      CALL message(TRIM(routine), 'Assuming a nested constellation in which an encompassed tile may sourround 1 tile itself')
    ELSE
      CALL message(TRIM(routine), 'Nested constellations are being ignored and geomtry is calculated exclusivly based on the characteristics of two tiles)')
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

    ! Create vertical axis
    tile_layer  => new_vgrid('connecttiles', ZAXIS_GENERIC, ncontis, levels=(/ (REAL(i,kind=wp),i=1,ncontis) /), units='')
    CALL register_vgrid(tile_layer)


  END SUBROUTINE Init_lex_config

#endif
END MODULE mo_lex_config_class
