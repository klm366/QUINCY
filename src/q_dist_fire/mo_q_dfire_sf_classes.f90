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

  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
  USE mo_q_dfire_sf_adaptor_class,  ONLY: sf_convert_units
  USE mo_q_dfire_sf_constants,  ONLY: df_sf_fc_1hr_idx, df_sf_fc_10hr_idx, df_sf_fc_100hr_idx, df_sf_fc_1000hr_idx !< fuel classes indices

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_sf_ros, t_sf_dynamics, t_sf_tree_damage, t_sf_ignitions, t_sf_mf, t_sf_lfmc, create_t_sf_prob_mort, &
            t_sf_atm, t_sf_fuel_load, t_sf_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_classes'

  !< Atmosphere related variables
  TYPE :: t_sf_atm
    !< To spitfire units (derived)
    REAL(wp), ALLOCATABLE :: t_air_min_c(:) ! [C]
    REAL(wp), ALLOCATABLE :: t_air_max_c(:) ! [C]
    REAL(wp), ALLOCATABLE :: wind_sf(:) ! [m/min]
    !< input
    REAL(wp), POINTER :: precip_daily_mmd(:) => NULL() ! [mm/day]
  END TYPE t_sf_atm

  INTERFACE t_sf_atm
    MODULE PROCEDURE :: create_sf_atm
  END INTERFACE t_sf_atm

  !< Rate of spread characteristics type
  TYPE t_sf_ros
    !< input
    INTEGER               :: nc
    REAL(wp), POINTER :: fbd(:) => NULL()                !< Weighted average across PFT for fuel bulk density [kg m-3]
    REAL(wp), POINTER :: fwd_wind_factor(:) => NULL()    !< Forward wind factor [-]

    !< outputs (intermediate)
    REAL(wp), ALLOCATABLE :: u_f_wind(:)                 !< Forward wind factor [-]
    REAL(wp), POINTER :: ros_forward(:) => NULL()        !< Forward rate of spread of the fire [m min-1]
    REAL(wp), POINTER :: ros_backwards(:) => NULL()      !< Backwards rate of spread of the fire [m min-1]
    REAL(wp), POINTER :: distance_travelled(:) => NULL() !< Total distance travelled [m]
    REAL(wp), POINTER :: l2b_ratio(:) => NULL()          !< length to breadth ratio [0 to 1]
    REAL(wp), POINTER :: fire_duration(:) => NULL()      !< [min]
    REAL(wp), POINTER :: mean_fire_area(:) => NULL()     !< [ha]

    REAL(wp), POINTER :: reaction_intensity(:) => NULL() !< I_R: the energy release rate per unit area of fire front [kJ m-2 min-1]
    REAL(wp), POINTER :: prop_flux_ratio(:) => NULL()    !< the proportion of IR that heats adjacent fuel particles to ignition [-]
    REAL(wp), POINTER :: wind_factor(:) => NULL()        !< wind_factor is a multiplier that accounts for the effect of wind in increasing the effective value of prop_flux_ratio[-]
    REAL(wp), POINTER :: eff_heating_number(:) => NULL() !< the proportion of a fuel particle that is heated to ...
                                                         !< the proportion of a fuel particle that is heated to [-]
    REAL(wp), POINTER :: pre_ign_heat(:) => NULL()       !< the amount of heat required to ignite a given mass of fuel [kJ kg-1]

    REAL(wp), POINTER :: opt_react_velo(:) => NULL()     !< Optimum reaction velocity [min-1]
    REAL(wp), POINTER :: net_fuel_load(:) => NULL()      !< Amount of fuel after substraction of its mineral content [kg m-2]
    REAL(wp), POINTER :: moist_damp_coeff(:) => NULL()   !< Moisture dampening coefficient [-]
    REAL(wp), POINTER :: mineral_damp_coeff(:) => NULL() !< Mineral dampening coefficient [-]
  END TYPE t_sf_ros

  INTERFACE t_sf_ros
    MODULE PROCEDURE :: create_sf_ros_init
  END INTERFACE t_sf_ros

  !<
  !< Calculated Spitfire fire dynamics
  !< Move spitfire dynamics data from calc_sf_dynamics <-> update_q_dfire_sf_dynamics subroutines
  !<
  TYPE t_sf_dynamics
    REAL(wp), POINTER :: fdi(:) => NULL()                     !< fire danger index [0 to 1]
    REAL(wp), POINTER :: area_burned(:) => NULL()             !< total area burned per unit area [ha ts-1]
    REAL(wp), POINTER :: area_burned_fract(:) => NULL()       !< fraction of area burned per unit area [0 to 1 ts-1]
    REAL(wp), POINTER :: dead_frac_consump(:,:) => NULL()     !< Dead fraction of consumption (#nc, #fuel classes) [0 to 1]
    REAL(wp), POINTER :: i_surface(:) => NULL()               !< Fire intensity [kW m-1]

    REAL(wp), ALLOCATABLE :: fc_dead_cls(:,:,:)               !< dead fuel consumption (#nc, #dead pools, #fuels classes) [g DM m-2]
    LOGICAL, ALLOCATABLE  :: is_fire_intense_enough(:)        !< mask, TRUE= fire intensity is strong enough to exist
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


  TYPE t_sf_ignitions
    ! inputs
    REAL(wp), POINTER :: a_nd(:) => NULL() ! [ignitions individuals-1 d-1]
    REAL(wp), POINTER :: pop_dens(:) => NULL() ! [individuals km-2]
    REAL(wp), POINTER :: lightning_freq(:) => NULL() ! [lightning ha-1 d-1]

    ! outputs
    REAL(wp), POINTER :: e_n_humans(:) => NULL() ! Rate of human caused ignitions [ha-1 ts-1]
    REAL(wp), POINTER :: e_n_lightning(:) => NULL() ! Rate of lightning caused ignitions [ha-1 ts-1]
    REAL(wp), POINTER :: e_n_ignitions(:) => NULL() ! Total rate of ignitions [ha-1 ts-1]

    ! config
    REAL(wp) :: max_human_pop_dens !< maximum human population density [individuals km-2]
  END TYPE t_sf_ignitions

  INTERFACE t_sf_ignitions
    MODULE PROCEDURE :: create_sf_ignitions
  END INTERFACE

  !> Livegrass fuel moisture content paremeters from Oberhagemann et al. (2025)
  TYPE :: t_sf_lfmc
    REAL(wp) :: min, max, b, c
  END TYPE t_sf_lfmc

  INTERFACE t_sf_lfmc
    MODULE PROCEDURE :: create_lfmc_from_lctlib
  END INTERFACE

  TYPE :: t_sf_mf !< fuel moisture
    REAL(wp), POINTER :: moisture_extinction(:) => NULL()
    ! REAL(wp), ALLOCATABLE :: val(:) !< TODO, omega_o (thonicke) and val (oberhagemann) shared

    !< Fuel moisture for Thonicke et al. (2010)
    INTEGER :: nalphas
    REAL(wp), ALLOCATABLE :: alphas(:) !< drying parameters
    REAL(wp), POINTER :: ab_moisture(:) => NULL() !< Above ground relative soil moisture layer [-]
    REAL(wp), POINTER :: omegas_thon(:,:) => NULL() !< 1 hr fuel moisture mixed with live grass(#nc)
    REAL(wp), ALLOCATABLE :: fuel_moisture_for_fc(:,:) !< fuel moisture used for fraction of dead fuel consumption (#nc, #dead fuel classes -1, 10, 100, 1000 hours-)
    !REAL(wp), ALLOCATABLE :: val(:) !< omega_o (mixed -1+lvg,10,100-) ! todo:

    !< Fuel moisture for oberhagemann et al. (2025)
    INTEGER :: sigma_dd_nr
    INTEGER :: sigma_dd_ros_nr
    INTEGER :: sigma_lv_nr
    REAL(wp), ALLOCATABLE :: alphas_dead(:) !<
    REAL(wp), ALLOCATABLE :: sigmas_dead(:) !<
    REAL(wp), ALLOCATABLE :: sigmas_live(:) !<
    REAL(wp), POINTER :: lvg_moisture(:) => NULL() !< livegrass fuel moisture [30% to 250%]
    REAL(wp), ALLOCATABLE :: dead(:,:) !< dead fuel moisture (#nc,#dead fuel moisture -1hr, 10hrs, 100hrs, 1000hrs, dry livegrass- )
    REAL(wp), ALLOCATABLE :: live(:,:) !< live fuel moisture (#nc, -wet livegrass, woody-)
    !< Weighted fuel moisture
    REAL(wp), ALLOCATABLE :: val(:) !< weighted fuel moisture (mixed-1,10,100 and lvg-, fine fuels-1hr and lvg-)
  END TYPE t_sf_mf

  INTERFACE t_sf_mf
    MODULE PROCEDURE :: create_new_sf_mf
  END INTERFACE


  !< Fuel load
  TYPE :: t_sf_fuel_load
    INTEGER :: nc
    INTEGER :: nfcl
    INTEGER :: nddpools, nlvpools
    REAL(wp), ALLOCATABLE :: dead_fuel_fcl(:,:,:) !< dead fuel (#nc, #dead fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), ALLOCATABLE :: cured_live_grass_fuel_fcl(:,:,:) !< cured live grass fuel (#nc, #dead fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), ALLOCATABLE :: live_grass_fuel_fcl(:,:,:) !< live grass fuel (#nc, #live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), ALLOCATABLE :: live_tree_fuel_fcl(:,:,:) !< live tree fuel (#nc, #live fuel pools, #fuel classes) [g DM m-2]
  END TYPE t_sf_fuel_load

  INTERFACE t_sf_fuel_load
    MODULE PROCEDURE :: create_sf_fuel_load
  END INTERFACE


  !< Configuration for the surface fire model
  TYPE :: t_sf_config
    INTEGER :: nc !< number of grid cells
    INTEGER :: nfcl_ros !< number of fuel classes involved in fire spread
    INTEGER :: nfcl !< total number of fuel classes
    INTEGER :: nddpools !< number of dead fuel pools

    REAL(wp) :: dtime !< time step length [s]
    REAL(wp) :: one_day !< length of one day in sec [s]

    INTEGER :: ros_model_id !< Select ROS implementation

    REAL(wp) :: const_calorific_heat_content !< Fuel calorific heat content [kJ kg-1]
    REAL(wp) :: const_isurface_extinct_thd !< Fire intensity threshold [kW m-1]
  END TYPE t_sf_config

  INTERFACE t_sf_config
    MODULE PROCEDURE :: create_sf_config
  END INTERFACE t_sf_config

  CONTAINS


  ! ======================================================================================================= !
  !>
  !> Initialize t_sf_ros
  !>
  FUNCTION create_sf_ros_init(nc, fwd_wind_factor, fuel_bulk_density, &
                             ros_forward, ros_backwards, &
                             distance_travelled, l2b_ratio, fire_duration, mean_fire_area, &
                             reaction_intensity, prop_flux_ratio, &
                             wind_factor, eff_heating_number, &
                             pre_ign_heat, &
                             opt_react_velo, net_fuel_load, &
                             moist_damp_coeff, mineral_damp_coeff) RESULT(sf_ros)
    INTEGER, INTENT(in)                  :: nc                      !< Total number of gridcells
    REAL(wp), POINTER, INTENT(in)        :: fwd_wind_factor(:)      !< Forward wind factor [-]
    REAL(wp), POINTER, INTENT(in)        :: fuel_bulk_density(:)    !< Weighted average across PFT for fuel bulk density [kg m-3]
    REAL(wp), POINTER, INTENT(inout)     :: ros_forward(:)          !<
    REAL(wp), POINTER, INTENT(inout)     :: ros_backwards(:)        !<
    REAL(wp), POINTER, INTENT(inout)     :: distance_travelled(:)   !<
    REAL(wp), POINTER, INTENT(inout)     :: l2b_ratio(:)            !<
    REAL(wp), POINTER, INTENT(inout)     :: fire_duration(:)        !<
    REAL(wp), POINTER, INTENT(inout)     :: mean_fire_area(:)       !<

    REAL(wp), POINTER, INTENT(inout)     :: reaction_intensity(:)   !<
    REAL(wp), POINTER, INTENT(inout)     :: prop_flux_ratio(:)      !<
    REAL(wp), POINTER, INTENT(inout)     :: wind_factor(:)          !<
    REAL(wp), POINTER, INTENT(inout)     :: eff_heating_number(:)   !<
    REAL(wp), POINTER, INTENT(inout)     :: pre_ign_heat(:)         !<
    REAL(wp), POINTER, INTENT(inout)     :: opt_react_velo(:)       !<
    REAL(wp), POINTER, INTENT(inout)     :: net_fuel_load(:)        !<
    REAL(wp), POINTER, INTENT(inout)     :: moist_damp_coeff(:)     !<
    REAL(wp), POINTER, INTENT(inout)     :: mineral_damp_coeff(:)   !<
    TYPE(t_sf_ros)                       :: sf_ros                  !< output [-]
    ! ----------------------------------------------------------------------------------------------------- !
    sf_ros%nc = nc
    sf_ros%fwd_wind_factor => fwd_wind_factor
    sf_ros%fbd => fuel_bulk_density ! input

    ! outputs
    sf_ros%ros_forward => ros_forward
    sf_ros%ros_backwards => ros_backwards
    sf_ros%distance_travelled => distance_travelled
    sf_ros%l2b_ratio => l2b_ratio
    sf_ros%fire_duration => fire_duration
    sf_ros%mean_fire_area => mean_fire_area

    sf_ros%reaction_intensity => reaction_intensity
    sf_ros%prop_flux_ratio => prop_flux_ratio
    sf_ros%wind_factor => wind_factor
    sf_ros%eff_heating_number => eff_heating_number
    sf_ros%pre_ign_heat => pre_ign_heat
    sf_ros%opt_react_velo => opt_react_velo
    sf_ros%net_fuel_load => net_fuel_load
    sf_ros%moist_damp_coeff => moist_damp_coeff
    sf_ros%mineral_damp_coeff => mineral_damp_coeff

    sf_ros%ros_forward(:) = 0._wp
    sf_ros%ros_backwards(:) = 0._wp
    sf_ros%distance_travelled(:) = 0._wp
    sf_ros%l2b_ratio(:) = 0._wp
    sf_ros%fire_duration(:) = 0._wp
    sf_ros%mean_fire_area = 0._wp

    sf_ros%reaction_intensity(:) = 0._wp
    sf_ros%prop_flux_ratio(:) = 0._wp
    sf_ros%wind_factor(:) = 0._wp
    sf_ros%eff_heating_number(:) = 0._wp
    sf_ros%pre_ign_heat(:) = 0._wp

    sf_ros%opt_react_velo(:) = 0._wp
    sf_ros%net_fuel_load(:) = 0._wp
    sf_ros%moist_damp_coeff(:) = 0._wp
    sf_ros%mineral_damp_coeff(:) = 0._wp

    ! derived
    ALLOCATE(sf_ros%u_f_wind(nc))
    sf_ros%u_f_wind(:) = 0._wp
  END FUNCTION create_sf_ros_init


  ! ======================================================================================================= !
  !>
  !> t_sf_atm constructor: Allocate and initialize
  !>
  FUNCTION create_sf_atm(t_air_min, t_air_max, wind_10m, precip_daily, precip_daily_mmd) RESULT (sf_atm)
    TYPE(t_sf_atm) :: sf_atm !< Output
    REAL(wp), INTENT(in) :: t_air_min(:) !< [K]
    REAL(wp), INTENT(in) :: t_air_max(:) !< [K]
    REAL(wp), INTENT(in) :: wind_10m(:) !< [m/min]
    REAL(wp), INTENT(in) :: precip_daily(:) !< in quincy units [kg m-2 s-1]
    REAL(wp), POINTER, INTENT(inout) :: precip_daily_mmd(:) !< in Spitfire units [mm/day]

    INTEGER :: nc

    sf_atm%precip_daily_mmd => precip_daily_mmd

    nc = SIZE(t_air_min)
    ALLOCATE(sf_atm%t_air_max_c(nc), sf_atm%t_air_min_c(nc), sf_atm%wind_sf(nc))

    sf_atm%t_air_max_c = sf_convert_units("Kelvin", "Celsius", t_air_max)
    sf_atm%t_air_min_c = sf_convert_units("Kelvin", "Celsius", t_air_min)
    sf_atm%precip_daily_mmd = sf_convert_units("kg m-2 s-1", "mm/day", precip_daily)
    sf_atm%wind_sf = sf_convert_units("m/s", "m/min", wind_10m) !< wind to spitfire units
  END FUNCTION create_sf_atm

  ! ======================================================================================================= !
  !>
  !> t_sf_dynamics constructor: Allocate and initialize
  !>
  FUNCTION create_sf_dynamics(nc, nfuels, ndeadpools, &
                              fdi, &
                              area_burned, area_burned_fract, &
                              dead_frac_consump, i_surface) &
                              RESULT (sf_dynamics)
    INTEGER, INTENT(in) :: nc          !< #Gridcells
    INTEGER, INTENT(in) :: nfuels      !< # fuel classes
    INTEGER, INTENT(in) :: ndeadpools  !< # dead fuel pools

    REAL(wp), POINTER, INTENT(in) :: fdi(:)                     !< fire danger index [0 to 1]
    REAL(wp), POINTER, INTENT(in) :: area_burned(:)             !< total area burned per unit area [ha ts-1]
    REAL(wp), POINTER, INTENT(in) :: area_burned_fract(:)       !< fraction of area burned per unit area [0 to 1 ts-1]
    REAL(wp), POINTER, INTENT(in) :: dead_frac_consump(:,:)     !< Dead fraction of consumption (#nc, #fuel classes) [0 to 1]
    REAL(wp), POINTER, INTENT(in) :: i_surface(:)               !< Fire intensity surface [kW m-1]

    TYPE(t_sf_dynamics) :: sf_dynamics !< Output

    sf_dynamics % fdi => fdi
    sf_dynamics % area_burned => area_burned
    sf_dynamics % area_burned_fract => area_burned_fract
    sf_dynamics % dead_frac_consump => dead_frac_consump
    sf_dynamics % i_surface => i_surface

    ALLOCATE(sf_dynamics % fc_dead_cls(nc, ndeadpools, nfuels)) ! spitfire dead fuel consumption (fuel classes)
    ALLOCATE(sf_dynamics % is_fire_intense_enough(nc))

    sf_dynamics % fdi = 0._wp
    sf_dynamics % area_burned = 0._wp
    sf_dynamics % area_burned_fract = 0._wp
    sf_dynamics % dead_frac_consump = 0._wp
    sf_dynamics % fc_dead_cls = 0._wp
    sf_dynamics % i_surface = 0._wp
    sf_dynamics % is_fire_intense_enough = .FALSE.

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


  ! ======================================================================================================= !
  !>
  !> Constructor for t_sf_ignitions: allocate and initialize
  !>
  FUNCTION create_sf_ignitions(nc, max_human_pop_dens, pop_dens, a_nd, lightning_freq, &
                        e_n_humans, e_n_lightning, e_n_ignitions) RESULT (sf_ignitions)
    INTEGER, INTENT(in)              :: nc
    REAL(wp), INTENT(in)             :: max_human_pop_dens
    REAL(wp), POINTER, INTENT(in)    :: pop_dens(:), a_nd(:), lightning_freq(:)
    REAL(wp), POINTER, INTENT(inout) :: e_n_humans(:), e_n_lightning(:), e_n_ignitions(:)
    TYPE(t_sf_ignitions)             :: sf_ignitions

    sf_ignitions%a_nd => a_nd
    sf_ignitions%pop_dens => pop_dens
    sf_ignitions%lightning_freq => lightning_freq
    sf_ignitions%max_human_pop_dens = max_human_pop_dens

    sf_ignitions%e_n_humans => e_n_humans
    sf_ignitions%e_n_lightning => e_n_lightning
    sf_ignitions%e_n_ignitions => e_n_ignitions

    sf_ignitions%e_n_humans = 0._wp
    sf_ignitions%e_n_lightning = 0._wp
    sf_ignitions%e_n_ignitions = 0._wp
  END FUNCTION create_sf_ignitions


  ! ======================================================================================================= !
  !>
  !> t_sf_mf fuel moisture constructor (Oberhagemann only)
  !>
  FUNCTION create_new_sf_mf(nc, nddfc, & ! core
                            moisture_extinction, &  ! core
                            nalphas_thon, alphas_thon, ab_moisture_thon, omegas_thon, &  ! thonicke only
                            sigma_dd_nr, sigma_dd_ros_nr, sigma_lv_nr, alphas_dead, sigmas_dead, sigmas_live, lvg_moisture &  ! oberhagemann only
                            ) RESULT (mf)
    INTEGER, INTENT(in) :: nc
    INTEGER, INTENT(in) :: nddfc !< number of dead fuel classes (1, 10, 100 and 1000 hours)
    REAL(wp), POINTER, INTENT(in) :: moisture_extinction(:)

    INTEGER, INTENT(in) :: nalphas_thon
    REAL(wp), INTENT(in) :: alphas_thon(:)
    REAL(wp), POINTER, INTENT(in) :: ab_moisture_thon(:)
    REAL(wp), POINTER, INTENT(inout) :: omegas_thon(:,:) !< omega_o and omega_nl, Thonicke only

    INTEGER, INTENT(in) :: sigma_dd_nr, sigma_dd_ros_nr, sigma_lv_nr
    REAL(wp), INTENT(in) :: alphas_dead(:)
    REAL(wp), INTENT(in) :: sigmas_dead(:), sigmas_live(:) !< Number of dead and live pools, Oberhagemann only
    REAL(wp), POINTER, INTENT(in) :: lvg_moisture(:)

    TYPE(t_sf_mf) :: mf

    ! core
    mf%moisture_extinction => moisture_extinction
    ALLOCATE(mf%fuel_moisture_for_fc(nc, nddfc))

    ! Thonicke et al. (2010)
    mf%nalphas = nalphas_thon
    ALLOCATE(mf%alphas, source = alphas_thon)
    mf%ab_moisture => ab_moisture_thon
    mf%omegas_thon => omegas_thon

    ! Oberhagemann et al. (2025)
    mf%sigma_dd_nr = sigma_dd_nr
    mf%sigma_lv_nr = sigma_lv_nr
    mf%sigma_dd_ros_nr = sigma_dd_ros_nr

    ALLOCATE(mf%alphas_dead, source = alphas_dead)
    ALLOCATE(mf%sigmas_dead, source = sigmas_dead)
    ALLOCATE(mf%sigmas_live, source = sigmas_live)
    mf%lvg_moisture => lvg_moisture

    ALLOCATE(mf%dead(nc, sigma_dd_nr))
    ALLOCATE(mf%live(nc, sigma_lv_nr))
    ALLOCATE(mf%val(nc))
    mf%dead = 0._wp
    mf%live = 0._wp
    mf%val = 0._wp
  END FUNCTION create_new_sf_mf


  ! ======================================================================================================= !
  !>
  !> t_sf_lfmc live grass moisture content constructor (Oberhagemann only)
  !>
  FUNCTION create_lfmc_from_lctlib(lctlib) RESULT (lfmc)
    TYPE(t_lctlib_element) :: lctlib
    TYPE(t_sf_lfmc) :: lfmc

    lfmc % min = lctlib%q_lfmc_min
    lfmc % max = lctlib%q_lfmc_max
    lfmc % b = lctlib%q_lfmc_b
    lfmc % c = lctlib%q_lfmc_c
  END FUNCTION create_lfmc_from_lctlib


  ! ======================================================================================================= !
  !>
  !> t_sf_fuel_load: initialize and rebuild fuel load
  !>
  FUNCTION create_sf_fuel_load(nc, ndd, nlv, nfcl, &
                              dead_fuel_1hr, dead_fuel_10hr, dead_fuel_100hr, dead_fuel_1000hr, &
                              live_tree_fuel_1hr, live_tree_fuel_10hr, live_tree_fuel_100hr, live_tree_fuel_1000hr, &
                              live_grass_fuel_1hr, cured_live_grass_fuel_1hr) RESULT (sf_fuel_load)
    INTEGER, INTENT(in) :: nc
    INTEGER, INTENT(in) :: ndd, nlv !< number of dead/live pools
    INTEGER, INTENT(in) :: nfcl !< number of fuel classes

    REAL(wp), INTENT(in) :: dead_fuel_1hr(:,:), dead_fuel_10hr(:,:), dead_fuel_100hr(:,:), dead_fuel_1000hr(:,:) !< dead fuel (#nc, #dead fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), INTENT(in) :: live_tree_fuel_1hr(:,:), live_tree_fuel_10hr(:,:), live_tree_fuel_100hr(:,:), live_tree_fuel_1000hr(:,:) !< live tree fuel (#nc, #live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), INTENT(in) :: live_grass_fuel_1hr(:,:) !< live grass fuel (#nc, #live fuel pools, #fuel classes) [g DM m-2]
    REAL(wp), INTENT(in) :: cured_live_grass_fuel_1hr(:,:) !< cured live grass fuel (#nc, #dead fuel pools, #fuel classes) [g DM m-2]

    TYPE(t_sf_fuel_load) :: sf_fuel_load

    sf_fuel_load%nc = nc
    sf_fuel_load%nddpools = ndd
    sf_fuel_load%nlvpools = nlv
    sf_fuel_load%nfcl = nfcl

    ALLOCATE(sf_fuel_load%dead_fuel_fcl(nc, ndd, nfcl))
    ALLOCATE(sf_fuel_load%cured_live_grass_fuel_fcl(nc, ndd, nfcl))
    ALLOCATE(sf_fuel_load%live_grass_fuel_fcl(nc, nlv, nfcl))
    ALLOCATE(sf_fuel_load%live_tree_fuel_fcl(nc, nlv, nfcl))

    sf_fuel_load%cured_live_grass_fuel_fcl = 0.0_wp
    sf_fuel_load%live_grass_fuel_fcl = 0.0_wp

    ! rebuild fuel variables
    sf_fuel_load%dead_fuel_fcl(:,:,df_sf_fc_1hr_idx) = dead_fuel_1hr ! dead fuel
    sf_fuel_load%dead_fuel_fcl(:,:,df_sf_fc_10hr_idx) = dead_fuel_10hr
    sf_fuel_load%dead_fuel_fcl(:,:,df_sf_fc_100hr_idx) = dead_fuel_100hr
    sf_fuel_load%dead_fuel_fcl(:,:,df_sf_fc_1000hr_idx) = dead_fuel_1000hr
    sf_fuel_load%live_grass_fuel_fcl(:,:,df_sf_fc_1hr_idx) = live_grass_fuel_1hr ! live grass fuel
    sf_fuel_load%cured_live_grass_fuel_fcl(:,:,df_sf_fc_1hr_idx) = cured_live_grass_fuel_1hr ! cured live grass fuel
    sf_fuel_load%live_tree_fuel_fcl(:,:,df_sf_fc_1hr_idx) = live_tree_fuel_1hr ! live tree fuel
    sf_fuel_load%live_tree_fuel_fcl(:,:,df_sf_fc_10hr_idx) = live_tree_fuel_10hr
    sf_fuel_load%live_tree_fuel_fcl(:,:,df_sf_fc_100hr_idx) = live_tree_fuel_100hr
    sf_fuel_load%live_tree_fuel_fcl(:,:,df_sf_fc_1000hr_idx) = live_tree_fuel_1000hr
  END FUNCTION create_sf_fuel_load


  ! ======================================================================================================= !
  !>
  !> t_sf_fuel_load: initialize and rebuild fuel load
  !>
  FUNCTION create_sf_config(nc, nfcl_ros, nfcl, nddpools, &
                            dtime, one_day, &
                            ros_model_id, &
                            const_calorific_heat_content, const_isurface_extinct_thd) &
                            RESULT(sf_config)
    INTEGER, INTENT(in) :: nc
    INTEGER, INTENT(in) :: nfcl_ros
    INTEGER, INTENT(in) :: nfcl
    INTEGER, INTENT(in) :: nddpools
    REAL(wp), INTENT(in) :: dtime
    REAL(wp), INTENT(in) :: one_day
    INTEGER, INTENT(in) :: ros_model_id
    REAL(wp), INTENT(in) :: const_calorific_heat_content
    REAL(wp), INTENT(in) :: const_isurface_extinct_thd
    TYPE(t_sf_config) :: sf_config

    sf_config % nc = nc
    sf_config % nfcl_ros = nfcl_ros
    sf_config % nfcl = nfcl
    sf_config % nddpools = nddpools
    sf_config % dtime = dtime
    sf_config % one_day = one_day
    sf_config % ros_model_id = ros_model_id
    sf_config % const_calorific_heat_content = const_calorific_heat_content
    sf_config % const_isurface_extinct_thd = const_isurface_extinct_thd
  END FUNCTION create_sf_config

#endif
END MODULE mo_q_dfire_sf_classes
