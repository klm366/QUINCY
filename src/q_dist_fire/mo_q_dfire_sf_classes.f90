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
!>#### Derived types for the fire model
!>
MODULE mo_q_dfire_sf_classes
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_sf_ros, t_sf_dynamics, t_sf_tree_damage

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_classes'

  !< Rate of spread characteristics type
  TYPE t_sf_ros
    INTEGER               :: nc
    REAL(wp), ALLOCATABLE :: ros_forward(:)        !< Forward rate of spread of the fire [m min-1]
    REAL(wp), ALLOCATABLE :: ros_backwards(:)      !< Backwards rate of spread of the fire [m min-1]
    REAL(wp), ALLOCATABLE :: distance_travelled(:) !< Total distance travelled [m]
    REAL(wp), ALLOCATABLE :: l2b_ratio(:)          !< length to breadth ratio [0 to 1]
    REAL(wp), ALLOCATABLE :: fire_duration(:)      !< [min]
    REAL(wp), ALLOCATABLE :: mean_fire_area(:)     !< [ha]

    REAL(wp), ALLOCATABLE :: reaction_intensity(:) !< I_R: the energy release rate per unit area of fire front [kJ m-2 min-1]
    REAL(wp), ALLOCATABLE :: prop_flux_ratio(:)    !< the proportion of IR that heats adjacent fuel particles to ignition [-]
    REAL(wp), ALLOCATABLE :: wind_factor(:)        !< wind_factor is a multiplier that accounts for the effect of wind in increasing the effective value of prop_flux_ratio[-]
    REAL(wp), ALLOCATABLE :: fbd(:)                !< Weighted average across PFT for fuel bulk density [kg m-3]
    REAL(wp), ALLOCATABLE :: eff_heating_number(:) !< the proportion of a fuel particle that is heated to ...
    !< the proportion of a fuel particle that is heated to [-]
    REAL(wp), ALLOCATABLE :: pre_ign_heat(:)       !< the amount of heat required to ignite a given mass of fuel [kJ kg-1]

    REAL(wp), ALLOCATABLE :: opt_react_velo(:)     !< Optimum reaction velocity [min-1]
    REAL(wp), ALLOCATABLE :: net_fuel_load(:)      !< Amount of fuel after substraction of its mineral content [kg m-2]
    REAL(wp), ALLOCATABLE :: moist_damp_coeff(:)   !< Moisture dampening coefficient [-]
    REAL(wp), ALLOCATABLE :: mineral_damp_coeff(:) !< Mineral dampening coefficient [-]
  END TYPE t_sf_ros

  INTERFACE t_sf_ros
    MODULE PROCEDURE :: create_sf_ros_init
  END INTERFACE t_sf_ros

  !<
  !< Calculated Spitfire fire dynamics
  !< Move spitfire dynamics data from calc_sf_dynamics <-> update_q_dfire_sf_dynamics subroutines
  !<
  TYPE t_sf_dynamics
    REAL(wp), ALLOCATABLE :: e_n_ignitions(:)
    REAL(wp), ALLOCATABLE :: e_n_humans(:)
    REAL(wp), ALLOCATABLE :: e_n_lightning(:)

    REAL(wp), ALLOCATABLE :: fdi(:) !< fire danger index [0 to 1]
    TYPE(t_sf_ros)        :: ros    !< fire dynamics

    REAL(wp), ALLOCATABLE :: area_burned(:)             !<
    REAL(wp), ALLOCATABLE :: area_burned_fract(:)

    REAL(wp), ALLOCATABLE :: dead_frac_consump(:,:)     !< Dead fraction of consumption (#nc, #fuel classes) [0 to 1]
    REAL(wp), ALLOCATABLE :: fc_dead_cls(:,:,:)         !< dead fuel consumption (#nc, #dead pools, #fuels classes) [g DM m-2]
    REAL(wp), ALLOCATABLE :: i_surface(:)               !< Fire intensity [kW m-1]

    LOGICAL, ALLOCATABLE  :: is_fire_intense_enough(:)  !< mask, TRUE= fire intensity is strong enough to exist
  END TYPE t_sf_dynamics

  INTERFACE t_sf_dynamics
    MODULE PROCEDURE :: create_sf_dynamics
  END INTERFACE t_sf_dynamics


  TYPE t_sf_tree_damage
    REAL(wp), ALLOCATABLE :: pm_ck(:)         !< probability of post fire mortality due to crown damage [-]
    REAL(wp), ALLOCATABLE :: r_ck             !< Resistance factor against crown damage [-]
    REAL(wp), ALLOCATABLE :: ck_p             !< Parameter for woody PFTs used in Pm (CK) equation [-]
    REAL(wp), ALLOCATABLE :: F_scorch_height  !< PFT-parameter in crown scorch equation [-]
    REAL(wp), POINTER :: pm_fire(:) => NULL() !< mix of tau and ck [-]
    REAL(wp), POINTER :: pm_tau(:)  => NULL() !< probability of post fire mortality due to cambial damage [-]
    REAL(wp), POINTER :: tau_l(:)             !< Residence time of the fire [min]
    REAL(wp), POINTER :: tau_c(:)             !< Critical time for cambial damage [min]
    REAL(wp), POINTER :: bark_thickness(:)    !< Tree barth thickness [cm]
    REAL(wp), POINTER :: prop_ck(:)           !< Fraction of the crown scorch [-]
    REAL(wp), POINTER :: scorch_height(:)     !< Scorch height [m]
    REAL(wp), POINTER :: crown_length(:)      !< Crown length of woody PFT [m]
  CONTAINS
    PROCEDURE :: copy_from_prob_mort
    GENERIC, PUBLIC :: assignment(=) => copy_from_prob_mort
  END TYPE t_sf_tree_damage


  INTERFACE t_sf_tree_damage
    MODULE PROCEDURE :: create_t_sf_prob_mort
  END INTERFACE t_sf_tree_damage
CONTAINS


  ! ======================================================================================================= !
  !>
  !> Initialize t_sf_ros
  !>
  FUNCTION create_sf_ros_init(nc) RESULT(sf_ros)
    INTEGER, INTENT(in)                  :: nc        !< Total number of gridcells
    TYPE(t_sf_ros)                       :: sf_ros    !< output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    sf_ros%nc = nc
    ALLOCATE(sf_ros%ros_forward(nc))
    ALLOCATE(sf_ros%ros_backwards(nc))
    ALLOCATE(sf_ros%distance_travelled(nc))
    ALLOCATE(sf_ros%l2b_ratio(nc))
    ALLOCATE(sf_ros%fire_duration(nc))
    ALLOCATE(sf_ros%mean_fire_area(nc))

    ALLOCATE(sf_ros%reaction_intensity(nc))
    ALLOCATE(sf_ros%prop_flux_ratio(nc))
    ALLOCATE(sf_ros%wind_factor(nc))
    ALLOCATE(sf_ros%fbd(nc))
    ALLOCATE(sf_ros%eff_heating_number(nc))
    ALLOCATE(sf_ros%pre_ign_heat(nc))

    ALLOCATE(sf_ros%opt_react_velo(nc))
    ALLOCATE(sf_ros%net_fuel_load(nc))
    ALLOCATE(sf_ros%moist_damp_coeff(nc))
    ALLOCATE(sf_ros%mineral_damp_coeff(nc))

    sf_ros%ros_forward(:) = 0._wp
    sf_ros%ros_backwards(:) = 0._wp
    sf_ros%distance_travelled(:) = 0._wp
    sf_ros%l2b_ratio(:) = 0._wp
    sf_ros%fire_duration(:) = 0._wp
    sf_ros%mean_fire_area = 0._wp

    sf_ros%reaction_intensity(:) = 0._wp
    sf_ros%prop_flux_ratio(:) = 0._wp
    sf_ros%wind_factor(:) = 0._wp
    sf_ros%fbd(:) = 0._wp
    sf_ros%eff_heating_number(:) = 0._wp
    sf_ros%pre_ign_heat(:) = 0._wp

    sf_ros%opt_react_velo(:) = 0._wp
    sf_ros%net_fuel_load(:) = 0._wp
    sf_ros%moist_damp_coeff(:) = 0._wp
    sf_ros%mineral_damp_coeff(:) = 0._wp
  END FUNCTION create_sf_ros_init


  ! ======================================================================================================= !
  !>
  !> t_sf_dynamics constructor: Allocate and initialize
  !>
  FUNCTION create_sf_dynamics(nc, nfuels, ndeadpools) RESULT (sf_spitfire)
    INTEGER, INTENT(in) :: nc          !< #Gridcells
    INTEGER, INTENT(in) :: nfuels      !< # fuel classes
    INTEGER, INTENT(in) :: ndeadpools  !< # dead fuel pools
    TYPE(t_sf_dynamics) :: sf_spitfire !< Output

    ALLOCATE(sf_spitfire % e_n_ignitions(nc))
    ALLOCATE(sf_spitfire % e_n_humans(nc))
    ALLOCATE(sf_spitfire % e_n_lightning(nc))
    ALLOCATE(sf_spitfire % fdi(nc))
    ALLOCATE(sf_spitfire % area_burned(nc))
    ALLOCATE(sf_spitfire % area_burned_fract(nc))
    ALLOCATE(sf_spitfire % dead_frac_consump(nc, nfuels))
    ALLOCATE(sf_spitfire % fc_dead_cls(nc, ndeadpools, nfuels)) ! spitfire dead fuel consumption (fuel classes)
    ALLOCATE(sf_spitfire % i_surface(nc))
    ALLOCATE(sf_spitfire % is_fire_intense_enough(nc))

    sf_spitfire % e_n_ignitions = 0._wp
    sf_spitfire % e_n_humans = 0._wp
    sf_spitfire % e_n_lightning = 0._wp
    sf_spitfire % fdi = 0._wp
    sf_spitfire % ros = t_sf_ros(nc)
    sf_spitfire % area_burned = 0._wp
    sf_spitfire % area_burned_fract = 0._wp
    sf_spitfire % dead_frac_consump = 0._wp
    sf_spitfire % fc_dead_cls = 0._wp
    sf_spitfire % i_surface = 0._wp
    sf_spitfire % is_fire_intense_enough = .FALSE.

  END FUNCTION create_sf_dynamics


  ! ======================================================================================================= !
  !>
  !> Constructor for t_sf_prob_mort
  !>
  FUNCTION create_t_sf_prob_mort(nc, pm_fire, pm_tau, tau_l, tau_c, &
                  bark_thickness, prop_ck, scorch_height, &
                  crown_length) RESULT (sf_pm)
    INTEGER, INTENT(in) :: nc
    REAL(wp), INTENT(inout), POINTER :: pm_fire(:)
    REAL(wp), INTENT(inout), POINTER :: pm_tau(:)
    REAL(wp), INTENT(inout), POINTER :: tau_l(:)
    REAL(wp), INTENT(inout), POINTER :: tau_c(:)
    REAL(wp), INTENT(inout), POINTER :: bark_thickness(:)
    REAL(wp), INTENT(inout), POINTER :: prop_ck(:)
    REAL(wp), INTENT(inout), POINTER :: scorch_height(:)
    REAL(wp), INTENT(inout), POINTER :: crown_length(:)
    TYPE(t_sf_tree_damage) :: sf_pm

    ALLOCATE(sf_pm%pm_ck(nc))
    ALLOCATE(sf_pm%r_ck)
    ALLOCATE(sf_pm%ck_p)
    ALLOCATE(sf_pm%F_scorch_height)
    sf_pm%pm_fire => pm_fire
    sf_pm%pm_tau => pm_tau
    sf_pm%tau_l => tau_l
    sf_pm%tau_c => tau_c
    sf_pm%bark_thickness => bark_thickness
    sf_pm%prop_ck => prop_ck
    sf_pm%scorch_height => scorch_height
    sf_pm%crown_length => crown_length

    sf_pm%pm_fire = 0._wp
    sf_pm%pm_tau = 0._wp
    sf_pm%pm_ck = 0._wp
    sf_pm%tau_l = 0._wp
    sf_pm%tau_c = 0._wp
    sf_pm%bark_thickness = 0._wp
    sf_pm%prop_ck = 0._wp
    sf_pm%r_ck = 0._wp
    sf_pm%ck_p = 0._wp
    sf_pm%scorch_height = 0._wp
    sf_pm%F_scorch_height = 0._wp
    sf_pm%crown_length = 0._wp
  END FUNCTION create_t_sf_prob_mort


  ! ======================================================================================================= !
  !>
  !> Copy from "other" into "self"
  !>
  SUBROUTINE copy_from_prob_mort(self, other)
    CLASS(t_sf_tree_damage), INTENT(inout) :: self
    CLASS(t_sf_tree_damage), INTENT(in)    :: other

    self % pm_ck = other % pm_ck ! copy data
    self % r_ck = other % r_ck
    self % ck_p = other % ck_p
    self % F_scorch_height = other % F_scorch_height
    self % pm_fire => other % pm_fire ! copy pointer reference
    self % pm_tau => other % pm_tau
    self % tau_l => other % tau_l
    self % tau_c => other % tau_c
    self % bark_thickness => other % bark_thickness
    self % prop_ck => other % prop_ck
    self % scorch_height => other % scorch_height
    self % crown_length => other % crown_length
  END SUBROUTINE copy_from_prob_mort

#endif
END MODULE mo_q_dfire_sf_classes
