!>
!! @brief Contains the memory class for the vegetation processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_memory_class

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_veg_memory, t_veg_phys_state, & !, t_veg_ecophys_state
            t_veg_water_state

  TYPE t_veg_phys_state
    REAL(wp), POINTER ::    &
      & lai            (:,  :), &  !< Leaf area index
      & fract_fpc      (:,  :), &  !< Foliage projected cover (fraction of ground area occupied by the vertical projection of foliage)
      & fract_fpc_max  (:,  :), &  !< Maximum fraction of tile that can be covered by leaves (assuming infinite lai)
      & fract_forest   (:,  :), &  !< Forest fraction
      & canopy_cond    (:,  :), &  !< Canopy conductance (1/resistance)
      & skin_can_depth (:,  :), &  !< Depth of skin reservoir on canopy
      & root_depth     (:,  :), &  !< Rooting depth [m]
      & root_depth_l   (:,:,:), &  !< Rooted depth per soil layer (until rooting depth) [m]
      & alb_vis_can    (:,  :), &  !< Albedo of the canopy in the visible range
      & alb_nir_can    (:,  :)     !< Albedo of the canopy in the NIR range
  END TYPE t_veg_phys_state

  TYPE t_veg_water_state
    REAL(wp), POINTER ::   &
      & wfract_can  (:,:), &  !< Wet (skin reservoir) fraction of canopy
      & wsr_can     (:,:), &  !< Water content in skin reservoir on canopy [m]
      & sfract_can  (:,:), &  !< Snow fraction of canopy
      & wsn_can     (:,:), &  !< Water content of snow reservoir onn canopy [m water equivalent]
      & ws_root     (:,:), &  !< Water content in root zone of the soil [m]
      & ws_fc_root  (:,:), &  !< Water content at field capacity in root zone of the soil [m]
      & ws_pwp_root (:,:), &  !< Water content at permanent wilting point in root zone of the soil [m]
      water_stress  (:,:)     !< Water stress factor (1: no water stress, 0: infinite stress)
  END TYPE t_veg_water_state

  TYPE t_veg_memory
    TYPE(t_veg_phys_state),  POINTER :: phys_state  => NULL()
    TYPE(t_veg_water_state), POINTER :: water_state => NULL()
  END TYPE t_veg_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'

END MODULE mo_veg_memory_class
