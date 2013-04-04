!>
!! @brief Contains types for the atmosphere-land interface memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-27)
!!
MODULE mo_atmland_memory_class

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_atm2land_memory, t_land2atm_memory

  TYPE t_atm2land_memory
    REAL(wp), POINTER :: &
      t_air            (:,:), &
      q_air            (:,:), &
      rain             (:,:), &
      snow             (:,:), &
      wind_air         (:,:), &
      wind_10m         (:,:), &
      lwrad_srf_down   (:,:), &
      swrad_srf_down   (:,:), &
      press_srf        (:,:), &
      drag_srf         (:,:), &
      t_acoef          (:,:), &
      t_bcoef          (:,:), &
      q_acoef          (:,:), &
      q_bcoef          (:,:), &
      cos_zenith_angle (:,:)
  END type t_atm2land_memory

  TYPE t_land2atm_memory
    REAL(wp), POINTER ::   &
      t_srf       (:,:),   &
      t_rad_srf   (:,:),   &
      qsat_srf    (:,:),   &
      fact_q_air  (:,:),   &
      fact_qsat_srf (:,:), &
      zh_srf      (:,:),   &
      zm_srf      (:,:),   &
      alb_vis_dif (:,:),   &
      alb_nir_dif (:,:),   &
      alb_vis_dir (:,:),   &
      alb_nir_dir (:,:)
  END type t_land2atm_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_memory_class'

END MODULE mo_atmland_memory_class
