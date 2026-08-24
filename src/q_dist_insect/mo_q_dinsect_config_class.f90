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

    !> Development flags
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
    !LOGICAL                       :: dinsect_bb_from_namelist        !< Bark beetle on/off damage data from namelist [string]
    !LOGICAL                       :: dinsect_df_from_namelist        !< Defoliator on/off damage data from namelist [string]

    !< Define dinsect_XX_from_namelist = FALSE to read from netcdf
    !CHARACTER(len=filename_max)   :: dinsect_bb_filename, &  !< Bark beetle data filename [string]
    !                                 dinsect_df_filename     !< Defoliator data filename [string]

    !< Prescribed for an entire year (daily data). For dev purposes only.
    !< Defined as 366 due to leap years. It is ignored for non-leap years.
    REAL(wp)          :: dinsect_dev_bb_forcing_dist_frac(366) !< Prescribed for an entire year (daily data). For dev purposes only.
    REAL(wp)          :: dinsect_dev_df_forcing_dist_frac(366) !< Prescribed for an entire year (daily data). For dev purposes only.
                                                     !< Todo: properly read netcdf files with daily data


    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':Init_dinsect_config'
    ! ----------------------------------------------------------------------------------------------------- !
    NAMELIST /lnd_q_dinsect_nml/     &
      active,                        &
      flag_dinsect,                  &
      dinsect_online,                &
      dinsect_dev_bb_forcing_dist_frac, &
      dinsect_dev_df_forcing_dist_frac

    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting dinsect configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                       = .TRUE.
    flag_dinsect                 = .FALSE.
    dinsect_online               = .FALSE.

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

    config%dinsect_dev_bb_forcing_dist_frac = dinsect_dev_bb_forcing_dist_frac
    config%dinsect_dev_df_forcing_dist_frac = dinsect_dev_df_forcing_dist_frac

    CALL close_nml(nml_handler)
  END SUBROUTINE Init_dinsect_config

#endif
END MODULE mo_q_dinsect_config_class
