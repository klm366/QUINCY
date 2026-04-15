!> Initialize insect disturbance variables (QUINCY)
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2019-2024, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md   Yimian Ma for dinsect
!> See LICENSES/ for license information
!> SPDX-License-Identifier: GPL-3.0-only
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### declare and define insect disturbance constants
!>
MODULE mo_q_dinsect_constants

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! Parameters for bark beetle effects
  REAL(wp), SAVE :: &
      tau_dead_leaf      , & !< Decay time for bark beetle impacted dead leaf falling down to groud [year]
      tau_dead_sapwood   , & !< Decay time for bark beetle impacted dead sapwood falling down to groud [year]
      tau_dead_heartwood     !< Decay time for bark beetle impacted dead heartwood falling down to groud [year]

  ! Parameters for defoliator effects
  REAL(wp), SAVE :: &
      defoli_eat_efficiency    , & !< Fraction of ingested VS defoliated CNP
      defoli_assim_frac_C        , & !< Fraction of assimilated VS ingested C
      defoli_assim_frac_N        , & !< Fraction of assimilated VS ingested N
      defoli_assim_frac_P        , & !< Fraction of assimilated VS ingested P
      defoli_growth_resp_frac  , & !< Fraction of growing-respired VS assimilated C
      !  , & !< Fraction of mass gained VS ingested C
      defoli_maint_resp_frac   , & !< Fraction of maintenance respirated C VS ingested C
      defoli_CP_ratio          , & !< the C P mass ratio of larva
      defoli_NP_ratio          , & !< the N P mass ratio of larva
      tau_insect                   !< Defoliator life span [day]

  ! Others
  REAL(wp), SAVE :: &
      C_mole_mass    , & !< Carbon molar mass
      N_mole_mass    , & !< Nitrogen molar mass
      P_mole_mass        !< Phosphorus molar mass

  ! Default values for damage fractions, only used when dinsect_from_namelist is OFF
  REAL(wp), SAVE :: &
      dinsect_const_damage_frac !< Default damage fraction

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_dinsect_constants'

CONTAINS

  ! ======================================================================================================= !
  !> Initialize constants used in the dinsect process
  !>
  SUBROUTINE init_dinsect_constants
    IMPLICIT NONE
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_dinsect_constants'

    ! Parameters for bark beetle effects
    tau_dead_leaf      = 3.0_wp      !< Decay time for leaf (years)
    tau_dead_sapwood   = 10.0_wp     !< Decay time for heartwood (years)
    tau_dead_heartwood = 10.0_wp     !< Decay time for heartwood (years)

    ! Parameters for defoliator effects
    defoli_eat_efficiency = 0.9_wp     !< Fraction of ingested VS defoliated CNP
    defoli_assim_frac_C = 0.5_wp         !< Fraction of assimilated VS ingested C
    defoli_assim_frac_N = 1.0_wp         !< Fraction of assimilated VS ingested N
    defoli_assim_frac_P = 1.0_wp         !< Fraction of assimilated VS ingested P
    defoli_growth_resp_frac = 0.45_wp  !< Fraction of growing-respiration VS assimilated C
    !defoli_mass_gained_frac = 0.3_wp  !< Fraction of mass gained VS ingested C
    defoli_maint_resp_frac = 0.0005_wp !< Fraction of maintenance respirated C VS insect C pool


    defoli_CP_ratio = 105.0_wp      !< C:P stoichiometric ratios of defoliator  mass ratio
    defoli_NP_ratio = 20.0_wp       !< N:P stoichiometric ratios of defoliator  mass ratio

    ! Insect life cycle parameters
    tau_insect = 45.0_wp            !< Insect lifespan (days)

    ! Molar masses
    C_mole_mass = 12.0_wp
    N_mole_mass = 14.0_wp
    P_mole_mass = 31.0_wp

    ! Default values for damage fractions, only used when dinsect_from_namelist is OFF
    dinsect_const_damage_frac = 0.001_wp

  END SUBROUTINE init_dinsect_constants

END MODULE mo_q_dinsect_constants
