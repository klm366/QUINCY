!> Insect disturbance process config (QUINCY)
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### define dinsect config structure, read dinsect namelist and init configuration parameters
!>
MODULE mo_q_dinsect_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: tolower
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_dinsect_config

  !-----------------------------------------------------------------------------------------------------
  !> Configuration of the insect disturbance process, derived from t_jsb_config
  !>
  !> Currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_q_dinsect_config
    !> Process control flags
    LOGICAL           :: flag_dinsect                 !< Enable insect disturbance process
    LOGICAL           :: dinsect_online               !< Enable coupling to vegetation and soil

    LOGICAL :: dinsect_flag_bb_from_namelist         !< Flag to indicate if disturbance BB fraction is read from namelist [logical]
    LOGICAL :: dinsect_flag_df_from_namelist         !< Flag to indicate if disturbance DF fraction is read from namelist [logical]
    !< Format:
    !< dinsect_bb_filename_prefix_XXXX.nc, where XXXX is the year
    CHARACTER(len=filename_max) :: dinsect_bb_filename_prefix !< Prefix for insect bb input file [string]
    CHARACTER(len=filename_max) :: dinsect_df_filename_prefix !< Prefix for insect df input file [string]
    CHARACTER(len=:), ALLOCATABLE :: dinsect_bb_path !< Path for insect bb input file [string]
    CHARACTER(len=:), ALLOCATABLE :: dinsect_df_path !< Path for insect bb input file [string]

    !> Development flags (leap year is ignored when not required)
    REAL(wp)          :: dinsect_dev_bb_forcing_dist_frac(366) !< Daily Prescribed mortality/defoliation rate for an entire year [0-1]
    REAL(wp)          :: dinsect_dev_df_forcing_dist_frac(366) !< Daily Prescribed mortality/defoliation rate for an entire year [0-1]

   CONTAINS
     PROCEDURE :: Init => Init_dinsect_config
  END type t_q_dinsect_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Configuration routine of t_q_dinsect_config
  !>
  !> Currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_dinsect_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,              ONLY: Register_vgrid
    USE mo_jsb_io,                ONLY: ZAXIS_GENERIC

    IMPLICIT NONE
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_dinsect_config), INTENT(inout) :: config    !< config type for dinsect
    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                       :: active
    CHARACTER(len=filename_max)   :: ic_filename, &
                                     bc_filename
    LOGICAL                       :: flag_dinsect                 !< on/off disturbance insect process [boolean]
    LOGICAL                       :: dinsect_online               !< on/off coupling dinsect to veg and soil [boolean]

    LOGICAL                       :: dinsect_flag_bb_from_namelist !< on/off reading bark beetle data from namelist [boolean]
    LOGICAL                       :: dinsect_flag_df_from_namelist !< on/off reading defoliator data from namelist [boolean]

    !< Used when dinsect_flag_bb_from_namelist or dinsect_flag_df_from_namelist is FALSE
    CHARACTER(len=filename_max)   :: dinsect_bb_filename_prefix !< prefix for bark beetle data filename [string]
    CHARACTER(len=filename_max)   :: dinsect_df_filename_prefix !< prefix for defoliator data filename [string]
    CHARACTER(len=:), ALLOCATABLE :: dinsect_bb_path !< Path for insect bb input file [string], must have / at the end
    CHARACTER(len=:), ALLOCATABLE :: dinsect_df_path !< Path for insect df input file [string], must have / at the end


    !< Used when dinsect_flag_df_from_namelist = TRUE
    !< Prescribed for an entire year (daily data). Every year will use the same forcing data. For dev purposes only.
    !< Defined as 366 due to leap years. It is ignored for non-leap years.
    REAL(wp)          :: dinsect_dev_bb_forcing_dist_frac(366) !< Prescribed for an entire year (daily data). For dev purposes only.
    REAL(wp)          :: dinsect_dev_df_forcing_dist_frac(366) !< Prescribed for an entire year (daily data). For dev purposes only.

    TYPE(t_jsb_vgrid), POINTER  :: time_doy_vgrid

    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':Init_dinsect_config'
    ! ----------------------------------------------------------------------------------------------------- !
    NAMELIST /lnd_q_dinsect_nml/     &
      active,                        &
      flag_dinsect,                  &
      dinsect_online,                &
      dinsect_flag_bb_from_namelist, &
      dinsect_flag_df_from_namelist, &
      dinsect_bb_filename_prefix,    &
      dinsect_df_filename_prefix,    &
      dinsect_bb_path,               &
      dinsect_df_path,               &
      dinsect_dev_bb_forcing_dist_frac, &
      dinsect_dev_df_forcing_dist_frac

    ! ----------------------------------------------------------------------------------------------------- !
    ! New Vgrid
    INTEGER :: i_doy
    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting dinsect configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                       = .TRUE.
    flag_dinsect                 = .FALSE.
    dinsect_online               = .FALSE.

    dinsect_flag_bb_from_namelist = .FALSE.
    dinsect_flag_df_from_namelist = .FALSE.
    dinsect_bb_filename_prefix    = 'bb_dist_frac_daily_'
    dinsect_df_filename_prefix    = 'df_dist_frac_daily_'
    dinsect_bb_path               = ''
    dinsect_df_path               = ''

    dinsect_dev_bb_forcing_dist_frac = 0._wp
    dinsect_dev_df_forcing_dist_frac = 0._wp

    ! ----------------------------------------------------------------------------------------------------- !
    ! read the namelist
    ! &lnd_q_dinsect_nml
    nml_handler = open_nml(TRIM(config%namelist_filename))
    ! ----------------------------------------------------------------------------------------------------- !
    nml_unit = position_nml('lnd_q_dinsect_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_q_dinsect_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! pass values as read from file
    config%active                       = active
    config%flag_dinsect                 = flag_dinsect
    config%dinsect_online               = dinsect_online

    config%dinsect_flag_bb_from_namelist = dinsect_flag_bb_from_namelist
    config%dinsect_flag_df_from_namelist = dinsect_flag_df_from_namelist
    config%dinsect_bb_filename_prefix    = dinsect_bb_filename_prefix
    config%dinsect_df_filename_prefix    = dinsect_df_filename_prefix
    config%dinsect_bb_path               = dinsect_bb_path
    config%dinsect_df_path               = dinsect_df_path

    config%dinsect_dev_bb_forcing_dist_frac = dinsect_dev_bb_forcing_dist_frac
    config%dinsect_dev_df_forcing_dist_frac = dinsect_dev_df_forcing_dist_frac

    ! create vertical grid
    time_doy_vgrid => new_vgrid('days_of_year', ZAXIS_GENERIC, 366, &
      & levels=(/ (REAL(i_doy, kind=wp), i_doy = 1,366) /), units='')
    CALL register_vgrid(time_doy_vgrid)

    CALL close_nml(nml_handler)
  END SUBROUTINE Init_dinsect_config

#endif
END MODULE mo_q_dinsect_config_class
