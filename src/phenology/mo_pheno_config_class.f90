!> Contains structures and methods for phenology config
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
MODULE mo_pheno_config_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pheno_config

  ! ======================================================================================================= !
  !> PHENO_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_pheno_config
    LOGICAL             :: l_forestRegrowth           !< If TRUE, the maxLAI is derived from biomass via allometric relationships
    CHARACTER(len=15)   :: scheme                     !< Phenology scheme to use ('logrop', 'climatology', 'external_lai')
    !LOGICAL :: first_call_phenology_this_day         !< TRUE marks if a new day has just started and
    !                                                 !< update_phenology is called for the top tile
    ! --- external_lai scheme settings ---
    CHARACTER(len=:), ALLOCATABLE :: extlai_mode             !< Temporal mode for external LAI: 'monthly' or 'annual'
    CHARACTER(len=:), ALLOCATABLE :: extlai_filename_prefix  !< Filename prefix for external LAI yearly files (prefix + YYYY.nc)
    REAL(wp)                      :: extlai_lai_max          !< Upper bound guard for external LAI values (m2/m2)
    CONTAINS
      PROCEDURE :: Init => Init_pheno_config
  END type t_pheno_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read PHENO_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_pheno_config(config)
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_pheno_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active, lrestart_cont
    LOGICAL                     :: l_forestRegrowth
    CHARACTER(len=15)           :: scheme
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    CHARACTER(len=10)           :: extlai_mode
    CHARACTER(len=filename_max) :: extlai_filename_prefix
    REAL(wp)                    :: extlai_lai_max
    NAMELIST /jsb_pheno_nml/  &
      & active,               &
      & lrestart_cont,        &
      & l_forestRegrowth,     &
      & scheme,               &
      & ic_filename,          &
      & bc_filename,          &
      & extlai_mode,          &
      & extlai_filename_prefix, &
      & extlai_lai_max
    INTEGER :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pheno_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting phenology configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active              = .TRUE.
    lrestart_cont       = .FALSE.      ! TRUE: Continue although PHENO variables are missing in restart file
    l_forestRegrowth    = .FALSE.
    scheme              = 'logrop'
    bc_filename         = 'bc_land_phenology.nc'
    ic_filename         = 'ic_land_phenology.nc'
    !first_call_phenology_this_day = .FALSE.
    extlai_mode         = 'monthly'
    extlai_filename_prefix = 'fc_land_lai_'
    extlai_lai_max      = 20.0_wp
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('jsb_pheno_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_pheno_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                     = active
    config%lrestart_cont              = lrestart_cont
    config%l_forestRegrowth           = l_forestRegrowth
    config%scheme                     = TRIM(scheme)
    config%ic_filename                = ic_filename
    config%bc_filename                = bc_filename
    config%extlai_mode                = TRIM(extlai_mode)
    config%extlai_filename_prefix     = TRIM(extlai_filename_prefix)
    config%extlai_lai_max             = extlai_lai_max

    CALL close_nml(nml_handler)

    ! Cross-check: external_lai scheme is not supported with QUINCY or with jsbach_forest_age_classes
    IF (TRIM(config%scheme) == 'external_lai') THEN
      IF (TRIM(config%model_config%model_scheme_char) == 'quincy') THEN
        CALL finish(TRIM(routine), &
          & 'The external_lai phenology scheme is not supported with QUINCY.')
      END IF
      IF (TRIM(config%model_config%usecase) == 'jsbach_forest_age_classes') THEN
        CALL finish(TRIM(routine), &
          & 'The external_lai phenology scheme is not supported with the jsbach_forest_age_classes usecase.')
      END IF
    END IF

  END SUBROUTINE Init_pheno_config

#endif
END MODULE mo_pheno_config_class
