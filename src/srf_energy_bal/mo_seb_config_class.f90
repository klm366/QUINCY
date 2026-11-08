!> Contains structures and methods for the configuration of the surface energy balance process
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
!>#### Configuration of the surface energy balance process
!>
!> The module defines the data structure [[t_seb_config]] that contains the configuration of the
!> surface energy balance (SEB_) process, i.e. all SEB variables that are configurable via namelist.
!> In subroutine [[Init_seb_config]] the namelist `jsb_seb_nml` is read and the configuration data
!> structure is set accordingly.
!>
MODULE mo_seb_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: lp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_config

  !> The data structure contains all configurable parameters of the SEB_ process. It extends
  !> [[mo_jsb_config_class:t_jsb_config]], a data structure containing basic parameters defined
  !> for all ICON-Land processes. Defaults are given in [[Init_seb_config]].
  TYPE, EXTENDS(t_jsb_config) :: t_seb_config
     CHARACTER(len=10)           :: scheme_albedo          !< Albedo scheme (default: 'echam5'), not used
     REAL(lp)                    :: lake_mixed_layer_depth !< Depth of lake mixed layer (default: 10) [m]
     REAL(lp)                    :: lake_min_ice_depth     !< Minimum ice thickness to start ice formation on
                                                           !< lakes (default: 0.05) [m]
     REAL(lp)                    :: coef_ril_tm1           !< Weight of the Richardson number at time t-1
                                                           !< (default: 0.15) []
     REAL(lp)                    :: coef_ril_t             !< Weight of the Richardson number at time t
                                                           !< (default: 0.70) []
     REAL(lp)                    :: coef_ril_tp1           !< Weight of the Richardson number at time t+1
                                                           !< (default: 0.15) []
                                                           !< The three weights are used in the standalone
                                                           !< model ([[mo_phy_schemes:update_drag]]) to
                                                           !< compute a drag coefficient that approximates the
                                                           !< one at time t but helps to maintain numerical
                                                           !< stability.
     REAL(lp)                    :: eb_threshold           !< Threshold for the residual of the surface energy
                                                           !< fluxes above which diagnostic messages are
                                                           !< printed by
                                                           !< [[mo_seb_interface:update_energy_balance]];
                                                           !< negative value: no printing
                                                           !< (default: -1) [W m-2]
     LOGICAL                     :: l_ice_on_lakes         !< Whether to use ice on lakes (default: true)
     INTEGER                     :: niter_tmx              !< For TMX: number of additional iterations of the
                                                           !< SEB solver, with exchange coefficients updated
                                                           !< in between (default: 0)
     LOGICAL                     :: l_skin_temp            !< Activate skin-temperature formulation (default:
                                                           !< false; not yet working with TMX)
   CONTAINS
     PROCEDURE :: Init => Init_seb_config
  END type t_seb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_config_class'

CONTAINS

  !> Read the namelist `jsb_seb_nml` and set the SEB_ configuration [[t_seb_config]].
  !>
  SUBROUTINE Init_seb_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_seb_config), INTENT(inout) :: config !< SEB_ configuration

    LOGICAL                     :: active, lrestart_cont            !< Namelist variables ([[t_seb_config]])
    CHARACTER(len=filename_max) :: ic_filename, bc_filename         !< Namelist variables
    CHARACTER(len=10)           :: scheme_albedo                    !< Namelist variable
    REAL(lp)                    :: &
      & lake_mixed_layer_depth, lake_min_ice_depth, eb_threshold, & !< Namelist variables
      & coef_ril_tm1, coef_ril_t, coef_ril_tp1                      !< Namelist variables
    LOGICAL                     :: l_ice_on_lakes                   !< Namelist variable
    INTEGER                     :: niter_tmx                        !< Namelist variable
    LOGICAL                     :: l_skin_temp                      !< Namelist variable

    NAMELIST /jsb_seb_nml/         &
      active,                      &
      lrestart_cont,               &
      ic_filename, bc_filename,    &
      scheme_albedo,               &
      lake_mixed_layer_depth,      &
      lake_min_ice_depth,          &
      coef_ril_tm1,                &
      coef_ril_t,                  &
      coef_ril_tp1,                &
      eb_threshold,                &
      l_ice_on_lakes,              &
      niter_tmx,                   &
      l_skin_temp

    INTEGER :: &
      & nml_handler, & !< Handler of the namelist file
      & nml_unit,    & !< Unit of the namelist file
      & istat          !< Status of positioning the namelist

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_seb_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting seb configuration')

    ! Set defaults
    active          = .TRUE.
    lrestart_cont   = .FALSE.  ! TRUE: Continue although SEB variables are missing in restart file
    ic_filename     = 'ic_land_srf.nc'
    bc_filename     = 'bc_land_srf.nc'
    scheme_albedo   = 'echam5'
    coef_ril_tm1    = 0.15_lp
    coef_ril_t      = 0.70_lp
    coef_ril_tp1    = 0.15_lp
    eb_threshold    = -1.0_lp  ! Negative values correspond to no printouts
    lake_mixed_layer_depth = 10._lp
    lake_min_ice_depth     = 0.05_lp
    l_ice_on_lakes  = .TRUE.
    niter_tmx       = 0
    l_skin_temp     = .FALSE.   ! default for QUINCY: TRUE

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_seb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_seb_nml)

    config%active                    = active
    config%lrestart_cont             = lrestart_cont
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%scheme_albedo             = scheme_albedo
    config%lake_mixed_layer_depth    = lake_mixed_layer_depth
    config%lake_min_ice_depth        = lake_min_ice_depth
    config%coef_ril_tm1              = coef_ril_tm1
    config%coef_ril_t                = coef_ril_t
    config%coef_ril_tp1              = coef_ril_tp1
    config%eb_threshold              = eb_threshold
    config%l_ice_on_lakes            = l_ice_on_lakes
    config%niter_tmx                 = niter_tmx
    config%l_skin_temp               = l_skin_temp

    CALL close_nml(nml_handler)

    IF (l_skin_temp) THEN
      CALL message(TRIM(routine), 'Skin temperature formulation is activated')
    END IF

  END SUBROUTINE Init_seb_config

#endif
END MODULE mo_seb_config_class
