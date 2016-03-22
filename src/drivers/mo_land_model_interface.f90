!!
!! @brief interface for land model  
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_land_model_interface

  USE mo_kind
  USE mo_base_memory

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: land_model_interface

  TYPE(TileClass), DIMENSION(ntile), SAVE        :: tile

  CONTAINS

  SUBROUTINE land_model_interface(t_air,q_air,rain,snow,wind_air,wind_10m,lw_srf_down,swvis_srf_down, &
                                  swnir_srf_down,swpar_srf_down,swpar_fdiffuse,press_srf,drag_srf,t_acoef,t_bcoef, & 
                                  q_acoef,q_bcoef,pch,cos_zenith_angle,co2_mixing_ratio, &
                                  ! undesirable fix for now
                                  veg_param) 

    USE mo_veg_parameters
    USE mo_veg_interface
    USE mo_veg_radiation
    USE mo_veg_photosynthesis

    IMPLICIT NONE
    ! 0.1 Input
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
      cos_zenith_angle , & 
      co2_mixing_ratio  
    TYPE(VegetationParameterClass) :: veg_param
    ! 0.4 Local
    REAL(wp) :: & 
      alb_vis_dif ,   &
      alb_nir_dif ,   &
      alb_vis_dir ,   &
      alb_nir_dir ,   &  
      ga, &
      vpd

      ga = 1._wp / ( drag_srf * wind_10m ) 

!write(*,*) 'interface C',tile%plant%pool%fine_root%carbon,&
! tile%plant%pool%fine_root%carbon,tile%plant%pool%coarse_root%carbon,&
! tile%plant%pool%sap_wood%carbon,tile%plant%pool%heart_wood%carbon,&
! tile%plant%pool%fruit%carbon,tile%plant%pool%labile%carbon
!write(*,*) 'interface N',tile%plant%pool%fine_root%nitrogen,& 
! tile%plant%pool%fine_root%nitrogen,tile%plant%pool%coarse_root%nitrogen,&
! tile%plant%pool%sap_wood%nitrogen,tile%plant%pool%heart_wood%nitrogen,&
! tile%plant%pool%fruit%nitrogen,tile%plant%pool%labile%nitrogen
!write(*,*) 'CN', tile%plant%pool%leaf%carbon/tile%plant%pool%leaf%nitrogen*12.0/14.0,&
!tile%plant%pool%fine_root%carbon/tile%plant%pool%fine_root%nitrogen*12.0/14.0 
tile(1)%plant%lai=6.2_wp
!tile(2)%plant%lai=5.0_wp
tile(1)%plant%pftid=1
!tile(2)%plant%pftid=2
     CALL update_radiation(tile)

     CALL update_canopy_fluxes(tile,veg_param)

     CALL update_vegetation(tile,t_air,veg_param)

  END SUBROUTINE land_model_interface

END MODULE mo_land_model_interface
