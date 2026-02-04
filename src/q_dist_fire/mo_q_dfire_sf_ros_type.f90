!> QUINCY
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
!> @brief Rate of spread derived type t_sf_ros
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_sf_ros_type
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_sf_ros

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_ros_type'

  !< Rate of spread characteristics type
  TYPE t_sf_ros
    INTEGER               :: nc
    REAL(wp), ALLOCATABLE :: ros_forward(:)        !< Forward rate of spread of the fire [m min-1]
    REAL(wp), ALLOCATABLE :: ros_backwards(:)      !< Backwards rate of spread of the fire [m min-1]
    REAL(wp), ALLOCATABLE :: distance_travelled(:) !< Total distance travelled [m]
    REAL(wp), ALLOCATABLE :: l2b_ratio(:)          !< length to breadth ratio [0 to 1]
    REAL(wp), ALLOCATABLE :: fire_duration(:)      !< [min]
    REAL(wp), ALLOCATABLE :: mean_fire_area(:)     !< [ha]

    REAL(wp), ALLOCATABLE :: reaction_intensity(:) !< I_R: the energy release rate per unit area of fire front [kJ m−2 min−1]
    REAL(wp), ALLOCATABLE :: prop_flux_ratio(:)    !< the proportion of IR that heats adjacent fuel particles to ignition [-]
    REAL(wp), ALLOCATABLE :: wind_factor(:)        !< wind_factor is a multiplier that accounts for the effect of wind in increasing the effective value of prop_flux_ratio[-]
    REAL(wp), ALLOCATABLE :: fbd(:)                !< Weighted average across PFT for fuel bulk density [kg m−3]
    REAL(wp), ALLOCATABLE :: eff_heating_number(:) !< the proportion of a fuel particle that is heated to ...
    !< the proportion of a fuel particle that is heated to [-]
    REAL(wp), ALLOCATABLE :: pre_ign_heat(:)       !< the amount of heat required to ignite a given mass of fuel [kJ kg−1]

    REAL(wp), ALLOCATABLE :: opt_react_velo(:)     !< Optimum reaction velocity [min-1]
    REAL(wp), ALLOCATABLE :: net_fuel_load(:)      !< Amount of fuel after substraction of its mineral content [kg m-2]
    REAL(wp), ALLOCATABLE :: moist_damp_coeff(:)   !< Moisture dampening coefficient [-]
    REAL(wp), ALLOCATABLE :: mineral_damp_coeff(:) !< Mineral dampening coefficient [-]
  END TYPE t_sf_ros

  INTERFACE t_sf_ros
    MODULE PROCEDURE :: create_sf_ros_init
  END INTERFACE t_sf_ros


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


#endif
END MODULE mo_q_dfire_sf_ros_type
