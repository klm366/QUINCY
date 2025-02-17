!> QUINCY simple harvest process config
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### define simple harvest config structure, read simple harvest namelist and init configuration parameters
!>
MODULE mo_q_sh_config_class
#ifndef __NO_QUINCY__

  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_math_constants,ONLY: eps4
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_q_sh_config

  ! ======================================================================================================= !
  !> Q_SH_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_q_sh_config
    LOGICAL   :: l_daily_harvest          !< determines if harvest is to be applied on a daily (TRUE) or annual (FALSE) basis
    LOGICAL   :: l_read_harvest_map       !< determines if harvest is to be read from a file (default FALSE), ignored in QS)
    REAL(wp)  :: harvest_fraction         !< default fraction used if not read from a map (annual change, relative to whole gridcell)
      !< for stand replacing harvest events, in QS, this is ignored and will be set to 1.0 instead

    ! QS specific configurations
    LOGICAL   :: harvest_active_in_qs     !< QS: if harvest should be active in QS (required for QS run script environment)
    LOGICAL   :: flag_stand_harvest       !< QS: if stand replacing harvest is enabled
    INTEGER   :: stand_replacing_year     !< QS: site specific stand replacing harvest year
    LOGICAL   :: flag_stand_harvest_event !< QS: if a site specific stand replacing harvest happened in this timestep
  CONTAINS
    PROCEDURE :: Init => Init_q_sh_config
  END TYPE t_q_sh_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sh_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read Q_SH_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_q_sh_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,              ONLY: Register_vgrid
    USE mo_jsb_io,                ONLY: ZAXIS_GENERIC
    USE mo_jsb_math_constants,    ONLY: eps4
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_sh_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename

    LOGICAL   :: l_daily_harvest          !< determines if harvest is to be applied on a daily (TRUE) or annual (FALSE) basis
    LOGICAL   :: l_read_harvest_map       !< determines if harvest is to be read from a file (default FALSE), ignored in QS)
    REAL(wp)  :: harvest_fraction         !< default fraction used if not read from a map (annual change, relative to whole gridcell)
      !< for stand replacing harvest events, in QS, this is ignored and will be set to 1.0 instead

    ! QS specific configurations
    LOGICAL   :: harvest_active_in_qs     !< QS: if harvest should be active in QS (required for QS run script environment)
    LOGICAL   :: flag_stand_harvest       !< QS: if stand replacing harvest is enabled
    INTEGER   :: stand_replacing_year     !< QS: site specific stand replacing harvest year

#ifdef __QUINCY_STANDALONE__
    NAMELIST /q_sh_ctl/  &
#else
    NAMELIST /q_sh_nml/  &
#endif
      & active,                   &
      & ic_filename,              &
      & bc_filename,              &
      & l_daily_harvest,          &
      & l_read_harvest_map,       &
      & harvest_fraction,         &
      & harvest_active_in_qs,     &
      & flag_stand_harvest,       &
      & stand_replacing_year


    INTEGER                     :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_sh_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting sh configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                        = .FALSE.
    bc_filename                   = 'bc_land_sh.nc'
    ic_filename                   = 'ic_land_sh.nc'
    !
    l_daily_harvest               = .TRUE.  ! determines if harvest is to be applied on a daily (TRUE) or annual (FALSE) basis
    l_read_harvest_map            = .FALSE. ! determines if harvest is to be read from a file (default FALSE, ignored in QS)
    ! for now arbitrary default fraction used if not read from a map (annual change, relative to whole gridcell)
    harvest_fraction              = 0.05_wp

    ! QS specific configurations
    harvest_active_in_qs          = .FALSE. ! (required for QS run script environment)
    flag_stand_harvest            = .FALSE.
    stand_replacing_year          = 1500
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
#ifdef __QUINCY_STANDALONE__
    nml_unit = position_nml('q_sh_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, q_sh_ctl)
#else
    nml_unit = position_nml('q_sh_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, q_sh_nml)
#endif

    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%l_daily_harvest           = l_daily_harvest
    config%l_read_harvest_map        = l_read_harvest_map
    config%harvest_fraction          = harvest_fraction
    config%harvest_active_in_qs      = harvest_active_in_qs
    config%flag_stand_harvest        = flag_stand_harvest
    config%stand_replacing_year      = stand_replacing_year

    ! QS only: Test if we have a stand replacing harvest event (default: FALSE) -- does not change in IQ!
    config%flag_stand_harvest_event  = .FALSE.

    CALL close_nml(nml_handler)

  END SUBROUTINE Init_q_sh_config

#endif
END MODULE mo_q_sh_config_class
