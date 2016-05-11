!>
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
  USE mo_exception,        ONLY: finish
  USE mo_base_memory,      ONLY: t_slm_tile
  USE mo_base_constants,   ONLY: ntile

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: land_model_interface

  TYPE(t_slm_tile), DIMENSION(ntile), SAVE        :: tile

  CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> land model interface
  !! 
  !! @todo land model interface is not working; science need to be integrated and improved
  !! @todo includes an undesirable fix for now (see subroutine input variables)
  !-----------------------------------------------------------------------------------------------------
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
    ! 0.1 InOut
    REAL(wp) ::          &
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
    TYPE(t_slm_veg_param) :: veg_param
    ! 0.2 Local
    REAL(wp) ::       & 
      alb_vis_dif ,   &
      alb_nir_dif ,   &
      alb_vis_dir ,   &
      alb_nir_dir ,   &  
      ga,             &
      vpd
    INTEGER  :: itile                        !< the current tile(i), for looping over ntile

!! stop here -----------------
CALL finish('land_model_interface','Stop! The land model is currently not working. Science needs to be improved to make it work.')
!! stop here -----------------


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

    DO itile = 1,ntile
      !tile%plant%pool%fine_root%carbon/tile%plant%pool%fine_root%nitrogen*12.0/14.0 
      tile(itile)%plant%lai=6.2_wp
      !tile(2)%plant%lai=5.0_wp
      tile(itile)%plant%pftid=1
      !tile(2)%plant%pftid=2

      CALL update_radiation(tile(itile))                         !< mo_veg_radiation
      CALL update_canopy_fluxes(tile(itile), veg_param)          !< mo_veg_photosynthesis
      CALL update_vegetation(tile(itile), veg_param)             !< mo_veg_interface
    END DO
  END SUBROUTINE land_model_interface

END MODULE mo_land_model_interface
