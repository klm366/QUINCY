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
    LOGICAL           :: dinsect_prescribe            !< Enable prescribed insect damage

    !> Input source flags
    LOGICAL           :: dinsect_from_namelist        !< Read damage data from namelist

    !> Insect parameters
    INTEGER           :: dinsect_type                 !< Insect type (1: Bark beetle, 2: Defoliator)
    INTEGER           :: dinsect_start_year           !< Start year for insect outbreak
    INTEGER           :: dinsect_year_span            !< Duration of outbreak [years]
    INTEGER           :: dinsect_start_doy            !< Start day of year for disturbance
    INTEGER           :: dinsect_duration             !< Duration within year [days]

    !> Impact parameters
    REAL(wp)          :: dinsect_damage_year1        !< Prescribed mortality/defoliation rate [0-1] for the first disturbance year
    REAL(wp)          :: dinsect_damage_year2        !< Prescribed mortality/defoliation rate [0-1] for the second disturbance year
    REAL(wp)          :: dinsect_damage_year3        !< Prescribed mortality/defoliation rate [0-1] for the third disturbance year
    REAL(wp)          :: dinsect_damage_year4        !< Prescribed mortality/defoliation rate [0-1] for the fourth disturbance year
    REAL(wp)          :: dinsect_damage_year5        !< Prescribed mortality/defoliation rate [0-1] for the fifth disturbance year
    REAL(wp)          :: dinsect_damage_year6        !< Prescribed mortality/defoliation rate [0-1] for the sixth disturbance year
    REAL(wp)          :: dinsect_damage_year7        !< Prescribed mortality/defoliation rate [0-1] for the seventh disturbance year
    REAL(wp)          :: dinsect_damage_year8        !< Prescribed mortality/defoliation rate [0-1] for the eighth disturbance year
    REAL(wp)          :: dinsect_damage_year9        !< Prescribed mortality/defoliation rate [0-1] for the ninth disturbance year
    REAL(wp)          :: dinsect_damage_year10       !< Prescribed mortality/defoliation rate [0-1]


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
    LOGICAL           :: flag_dinsect                 !< on/off disturbance insect process [boolean]
    LOGICAL           :: dinsect_online               !< on/off coupling dinsect to veg and soil [boolean]
    LOGICAL           :: dinsect_prescribe            !< on/off prescribing insect damage or not [boolean]
    LOGICAL           :: dinsect_from_namelist        !< on/off damage data from namelist [boolean]
    INTEGER           :: dinsect_type                 !< insect type. 1: Bark beetle [Default]; 2: Defoliator
    INTEGER           :: dinsect_start_year           !< prescribed start year for insect outbreak
    INTEGER           :: dinsect_year_span            !< the time span of insect outbreak (unit: years)
    INTEGER           :: dinsect_start_doy            !< the doy of starting day of a insect disturbance
    INTEGER           :: dinsect_duration             !< the time duration that disturbance happened during a year [unit: days]

    ! variables
    REAL(wp)          :: dinsect_damage_year1        !< Prescribed mortality/defoliation rate [0-1] for the first disturbance year
    REAL(wp)          :: dinsect_damage_year2        !< Prescribed mortality/defoliation rate [0-1] for the second disturbance year
    REAL(wp)          :: dinsect_damage_year3        !< Prescribed mortality/defoliation rate [0-1] for the third disturbance year
    REAL(wp)          :: dinsect_damage_year4        !< Prescribed mortality/defoliation rate [0-1] for the fourth disturbance year
    REAL(wp)          :: dinsect_damage_year5        !< Prescribed mortality/defoliation rate [0-1] for the fifth disturbance year
    REAL(wp)          :: dinsect_damage_year6        !< Prescribed mortality/defoliation rate [0-1] for the sixth disturbance year
    REAL(wp)          :: dinsect_damage_year7        !< Prescribed mortality/defoliation rate [0-1] for the seventh disturbance year
    REAL(wp)          :: dinsect_damage_year8        !< Prescribed mortality/defoliation rate [0-1] for the eighth disturbance year
    REAL(wp)          :: dinsect_damage_year9        !< Prescribed mortality/defoliation rate [0-1] for the ninth disturbance year
    REAL(wp)          :: dinsect_damage_year10       !< Prescribed mortality/defoliation rate [0-1] for the tenth disturbance year


    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':Init_dinsect_config'
    ! ----------------------------------------------------------------------------------------------------- !
    NAMELIST /lnd_q_dinsect_nml/       &
      flag_dinsect,                  &
      dinsect_online,                &
      dinsect_type,                  &
      dinsect_prescribe,             &
      dinsect_start_year,            &
      dinsect_year_span,             &
      dinsect_from_namelist,         &
      dinsect_duration,              &
      dinsect_start_doy,             &
      dinsect_damage_year1,          &
      dinsect_damage_year2,          &
      dinsect_damage_year3,          &
      dinsect_damage_year4,          &
      dinsect_damage_year5,          &
      dinsect_damage_year6,          &
      dinsect_damage_year7,          &
      dinsect_damage_year8,          &
      dinsect_damage_year9,          &
      dinsect_damage_year10

    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting dinsect configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    flag_dinsect                 = .FALSE.
    dinsect_online               = .FALSE.
    dinsect_type                 = 1         ! 1: Bark beetles 2: Defoliators
    dinsect_prescribe            = .FALSE.
    dinsect_from_namelist        = .TRUE.
    dinsect_start_year           = 2001
    dinsect_year_span            = 3
    dinsect_duration             = 180
    dinsect_start_doy            = 152

    dinsect_damage_year1          = 0._wp
    dinsect_damage_year2          = 0._wp
    dinsect_damage_year3          = 0._wp
    dinsect_damage_year4          = 0._wp
    dinsect_damage_year5          = 0._wp
    dinsect_damage_year6          = 0._wp
    dinsect_damage_year7          = 0._wp
    dinsect_damage_year8          = 0._wp
    dinsect_damage_year9          = 0._wp
    dinsect_damage_year10         = 0._wp

    ! ----------------------------------------------------------------------------------------------------- !
    ! read the namelist
    ! &lnd_q_dinsect_nml
    nml_handler = open_nml(TRIM(config%namelist_filename))
    ! ----------------------------------------------------------------------------------------------------- !
    nml_unit = position_nml('lnd_q_dinsect_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_q_dinsect_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! pass values as read from file
    config%flag_dinsect                 = flag_dinsect
    config%dinsect_online               = dinsect_online
    config%dinsect_type                 = dinsect_type
    config%dinsect_prescribe            = dinsect_prescribe
    config%dinsect_from_namelist        = dinsect_from_namelist
    config%dinsect_start_year           = dinsect_start_year
    config%dinsect_year_span            = dinsect_year_span
    config%dinsect_duration             = dinsect_duration
    config%dinsect_start_doy            = dinsect_start_doy
    config%dinsect_damage_year1         = dinsect_damage_year1
    config%dinsect_damage_year2         = dinsect_damage_year2
    config%dinsect_damage_year3         = dinsect_damage_year3
    config%dinsect_damage_year4         = dinsect_damage_year4
    config%dinsect_damage_year5         = dinsect_damage_year5
    config%dinsect_damage_year6         = dinsect_damage_year6
    config%dinsect_damage_year7         = dinsect_damage_year7
    config%dinsect_damage_year8         = dinsect_damage_year8
    config%dinsect_damage_year9         = dinsect_damage_year9
    config%dinsect_damage_year10        = dinsect_damage_year10

    CALL close_nml(nml_handler)
  END SUBROUTINE Init_dinsect_config

#endif
END MODULE mo_q_dinsect_config_class