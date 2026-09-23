! Author: Dr. Sri Nandini-Weiß
! Date: 2026-09-23
!> Coastal-constraint process configuration
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
!>#### Contains namelist information for the coastal-constraint process
!>
!> The process applies the salinity and tidal-inundation constraints on mangrove
!> productivity of Zhao et al. (2025, ORCHIDEE-MAN, doi:10.1029/2025MS005185).
!> Their equation numbers are used throughout this process.
!>
!> Two switches separate calculation from application, as for the trait process.
!> With active true and both application switches false the control factors are
!> computed and written as diagnostics but productivity is unchanged.
!>
MODULE mo_coast_config_class
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message
  USE mo_io_units,         ONLY: filename_max
  USE mo_jsb_config_class, ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_coast_config

  TYPE, EXTENDS(t_jsb_config) :: t_coast_config
    ! l_cs_one_sided added: sensitivity switch for the salinity factor.
    ! Author: Dr. Sri Nandini-Weiß
    ! Date: 2026-09-23
    LOGICAL  :: &
      & l_apply_salinity_stress, & !< Apply the salinity control factor to assimilation
      & l_apply_tidal_stress,    & !< Apply the tidal-inundation control factor to assimilation
      & l_normalise_cs,          & !< Scale the salinity factor so that its maximum is one
      & l_cs_one_sided             !< No salinity penalty below the optimum salinity
    REAL(wp) :: &
      & cs_min,                  & !< Minimum of the salinity control factor
      & cin_min,                 & !< Minimum of the tidal-inundation control factor
      & b_agrs_lim,              & !< Stilt-root biomass giving full ventilation [Mg ha-1]
      & b_agrp_lim,              & !< Pneumatophore biomass giving full ventilation [Mg ha-1]
      & vent_weight_s,           & !< Weight of the stilt-root ventilation term
      & vent_weight_p,           & !< Weight of the pneumatophore ventilation term
      & tree_density               !< Assumed stand density [trees ha-1]
  CONTAINS
    PROCEDURE :: Init => Init_coast_config
  END TYPE t_coast_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_coast_config_class'

CONTAINS

  ! ======================================================================================================== !
  !>
  !>#### Setup the coastal-constraint configuration
  !>
  !> We here configure the coastal process from namelist 'jsb_coast_nml'. The published
  !> parameter values of Zhao et al. (2025, their Table 1) are used as defaults and are
  !> exposed so that they can be changed without recompiling.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE Init_coast_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_coast_config), INTENT(inout) :: config

    LOGICAL  :: active, lrestart_cont
    ! l_cs_one_sided added to the local namelist variables.
    ! Author: Dr. Sri Nandini-Weiß
    ! Date: 2026-09-23
    LOGICAL  :: l_apply_salinity_stress, l_apply_tidal_stress, l_normalise_cs, l_cs_one_sided
    REAL(wp) :: cs_min, cin_min, b_agrs_lim, b_agrp_lim
    REAL(wp) :: vent_weight_s, vent_weight_p, tree_density
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    ! l_cs_one_sided added to jsb_coast_nml.
    ! Author: Dr. Sri Nandini-Weiß
    ! Date: 2026-09-23
    NAMELIST /jsb_coast_nml/     &
      & active,                  &
      & lrestart_cont,           &
      & ic_filename,             &
      & bc_filename,             &
      & l_apply_salinity_stress, &
      & l_apply_tidal_stress,    &
      & l_normalise_cs,          &
      & l_cs_one_sided,          &
      & cs_min,                  &
      & cin_min,                 &
      & b_agrs_lim,              &
      & b_agrp_lim,              &
      & vent_weight_s,           &
      & vent_weight_p,           &
      & tree_density

    INTEGER :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_coast_config'
    ! ------------------------------------------------------------------------------------------------------ !

    CALL message(TRIM(routine), 'Starting coastal configuration')

    active                  = .FALSE.
    lrestart_cont           = .FALSE.
    bc_filename             = 'bc_land_coastal.nc'
    ic_filename             = 'ic_land_coastal.nc'
    l_apply_salinity_stress = .FALSE.
    l_apply_tidal_stress    = .FALSE.
    ! The published quadratic reaches 0.8402 at its optimum rather than one. With
    ! l_normalise_cs false the equation is used exactly as published, so a penalty of
    ! about 16 percent is applied even at the most favourable salinity. Setting it true
    ! divides by that maximum and is a deliberate departure from the publication.
    l_normalise_cs          = .FALSE.
    ! Sensitivity switch: with l_cs_one_sided true, salinities below the optimum of the
    ! quadratic (29.907 PSS) are treated as the optimum, so only high salinity reduces
    ! productivity. The published low-salinity branch rests on few observations. This is
    ! a departure from the publication; the default keeps Equation 1 as published.
    ! Author: Dr. Sri Nandini-Weiß
    ! Date: 2026-09-23
    l_cs_one_sided          = .FALSE.
    cs_min                  = 0.5_wp     ! Zhao et al. (2025), Table 1
    cin_min                 = 0.5_wp     ! Zhao et al. (2025), Table 1
    b_agrs_lim              = 4.0_wp     ! Zhao et al. (2025), Section 2.2.3
    b_agrp_lim              = 0.8_wp     ! Zhao et al. (2025), Section 2.2.3
    vent_weight_s           = 0.5_wp     ! equal weighting, no global species map
    vent_weight_p           = 0.5_wp
    tree_density            = 1000.0_wp  ! trees ha-1

    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('jsb_coast_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_coast_nml)
    CALL close_nml(nml_handler)

    config%active                  = active
    config%lrestart_cont           = lrestart_cont
    config%ic_filename             = ic_filename
    config%bc_filename             = bc_filename
    config%l_apply_salinity_stress = l_apply_salinity_stress
    config%l_apply_tidal_stress    = l_apply_tidal_stress
    config%l_normalise_cs          = l_normalise_cs
    ! Store the one-sided salinity switch.
    ! Author: Dr. Sri Nandini-Weiß
    ! Date: 2026-09-23
    config%l_cs_one_sided          = l_cs_one_sided
    config%cs_min                  = cs_min
    config%cin_min                 = cin_min
    config%b_agrs_lim              = b_agrs_lim
    config%b_agrp_lim              = b_agrp_lim
    config%vent_weight_s           = vent_weight_s
    config%vent_weight_p           = vent_weight_p
    config%tree_density            = tree_density

    IF (active .AND. (l_apply_salinity_stress .OR. l_apply_tidal_stress)) &
      & CALL message(TRIM(routine), 'Coastal control factors will be applied to assimilation')

  END SUBROUTINE Init_coast_config

#endif
END MODULE mo_coast_config_class
