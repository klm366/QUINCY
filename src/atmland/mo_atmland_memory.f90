!!
!! @brief memory structure for atm->land and land->atm processes (not used yet) 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_atmland_memory

  USE mo_kind, ONLY: wp

  TYPE t_atm2land_memory

    REAL(wp) :: &
      t_air            , &
      q_air            , &
      rain             , &
      snow             , &
      wind_air         , &
      wind_10m         , &
      lw_srf_down      , &
      swvis_srf_down   , &
      swnir_srf_down   , &
      swpar_srf_down   , &
      swpar_fdiffuse   , &
      press_srf        , &
      drag_srf         , &
      t_acoef          , &
      t_bcoef          , &
      q_acoef          , &
      q_bcoef          , &
      pch              , &
      cos_zenith_angle

  END TYPE t_atm2land_memory

  TYPE t_land2atm_memory

    REAL(wp) :: &
      t_srf       ,   &
      t_eff_srf   ,   &
      qsat_srf    ,   &
      fact_q_air  ,   &
      fact_qsat_srf , &
      zh_srf      ,   &
      zm_srf      ,   &
      alb_vis_dif ,   &
      alb_nir_dif ,   &
      alb_vis_dir ,   &
      alb_nir_dir 

  END TYPE t_land2atm_memory

CONTAINS

END MODULE mo_atmland_memory
