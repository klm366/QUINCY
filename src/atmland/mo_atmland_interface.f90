!>
!! @brief Interface to run JSBACH for one time step
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Integrated into SLM by Jan Engel (20160823) using JSBACH4:git
!!
MODULE mo_atmland_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish
  USE mo_base_memory,             ONLY: t_slm_tile, t_jsb_task_options
  USE mo_jsb_var_class,           ONLY: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_real1d
  USE mo_jsb_process_class,       ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)   !  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  dsl4jsb_Use_memory(L2A_)   !  USE mo_l2a_memory_class, ONLY: t_l2a_memory

!!! jsbach4
!!  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
!!  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
!!  USE mo_jsb_config_class,   ONLY: t_jsb_config
!!  USE mo_jsb_process_class,  ONLY: t_jsb_process
!!  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
!!  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE, Contains_lct

!!  dsl4jsb_Use_processes A2L_, SEB_, TURB_, SOIL_, HYDRO_, RAD_
!!  dsl4jsb_Use_memory(A2L_)
!!  dsl4jsb_Use_memory(SEB_)
!!  dsl4jsb_Use_memory(TURB_)
!!  dsl4jsb_Use_memory(SOIL_)
!!  dsl4jsb_Use_memory(HYDRO_)
!!  dsl4jsb_Use_memory(RAD_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_a2l, update_l2a


  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> pass atmland variables from driver into a2l memory per timestep
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_a2l( &
    ! the sequence of input variables is consistent with the sequence of var in t_a2l_memory !
    ! some variables differ compared to the JSBACH4 routine, see comments 
    ! all variables are defined as: real(wp), optional, intent(in)
      tile               , &
      options            , &
      t_air              , &
      q_air              , &
      rain               , &
      snow               , &
      wind_air           , &
      wind_10m           , &
      lw_srf_down        , &
      swvis_srf_down     , &
      swnir_srf_down     , &
      swpar_srf_down     , &
      swpar_fdiffuse     , &  ! jsbach4: frac_par_diffuse
      press_srf          , &
      drag_srf           , &
      t_acoef            , &
      t_bcoef            , &
      q_acoef            , &
      q_bcoef            , &
      pch                , &
      cos_zenith_angle   , &
!      CO2_air            , &   ! jsbach4 variables, replaced by co2_mixing_ratio  _c13  _c14
!      CO2_air_mol        , &   !
      co2_mixing_ratio     , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
      co2_mixing_ratio_C13 , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
      co2_mixing_ratio_C14 , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
      ga                   , &  !> @todo ga in t_a2l_memory: is a variable from SLM, which may be moved elsewhere
      nhx_deposition     , &
      noy_deposition     , &
      p_deposition       , &
      ! For lakes:
      albedo_lwtr        , &
      drag_wtr           , &
      drag_ice           , &
      t_acoef_wtr        , &
      t_bcoef_wtr        , &
      q_acoef_wtr        , &
      q_bcoef_wtr        , &
      t_acoef_ice        , &
      t_bcoef_ice        , &
      q_acoef_ice        , &
      q_bcoef_ice          &
      & )

!    USE mo_jsb_physical_constants, ONLY: molarMassDryAir, molarMassCO2

    IMPLICIT NONE 
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< one tile
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    REAL(wp), OPTIONAL, INTENT(in) :: &
                                      t_air(:)              , &
                                      q_air(:)              , &
                                      rain(:)               , &
                                      snow(:)               , &
                                      wind_air(:)           , &
                                      wind_10m(:)           , &
                                      lw_srf_down(:)        , &
                                      swvis_srf_down(:)     , &
                                      swnir_srf_down(:)     , &
                                      swpar_srf_down(:)     , &
                                      swpar_fdiffuse(:)     , &  ! jsbach4: frac_par_diffuse
                                      press_srf(:)          , &
                                      drag_srf(:)           , &
                                      t_acoef(:)            , &
                                      t_bcoef(:)            , &
                                      q_acoef(:)            , &
                                      q_bcoef(:)            , &
                                      pch(:)                , &
                                      cos_zenith_angle(:)   , &
                                !      CO2_air(:)            , &   ! jsbach4 variables, replaced by co2_mixing_ratio  _c13  _c14
    ! Info: in JSBACH4 CO2_air_mol(:) is calculated as: dsl4jsb_var2D_onChunk(A2L_, CO2_air_mol) = CO2_air(:) * molarMassDryAir / molarMassCO2
                                      co2_mixing_ratio(:)     , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
                                      co2_mixing_ratio_C13(:) , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
                                      co2_mixing_ratio_C14(:) , &  ! from SLM, replacing CO2_air and CO2_air_mol from jsbach4
                                      ga(:)                   , &  !> @todo ga in t_a2l_memory: is a variable from SLM, which may be moved elsewhere
                                      nhx_deposition(:)     , &    ! added with r506
                                      noy_deposition(:)     , &    ! added with r506
                                      p_deposition(:)       , &    ! added with r... ticket #261
                                      ! For lakes:
                                      drag_wtr(:)           , &
                                      drag_ice(:)           , &
                                      t_acoef_wtr(:)        , &
                                      t_bcoef_wtr(:)        , &
                                      q_acoef_wtr(:)        , &
                                      q_bcoef_wtr(:)        , &
                                      t_acoef_ice(:)        , &
                                      t_bcoef_ice(:)        , &
                                      q_acoef_ice(:)        , &
                                      q_bcoef_ice(:)        , &
                                      albedo_lwtr(:)
    ! ---------------------------
    ! 0.2 Local
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_a2l'
!    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)   !     TYPE(t_a2l_memory), POINTER :: a2l__mem
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    IF (PRESENT(t_air)) THEN  ! not present for e.g. test_canopy model  *slmtemp
      IF (ice-ics+1 /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions') 
    END IF
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(A2L_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(A2L_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)       ! a2l__mem => a2l_memory




    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! 1.0 pass variables values to memory
    IF (PRESENT(t_air))             dsl4jsb_var2D_onChunk(A2L_,t_air)             = t_air             (:)
    IF (PRESENT(q_air))             dsl4jsb_var2D_onChunk(A2L_,q_air)             = q_air             (:)
    IF (PRESENT(rain))              dsl4jsb_var2D_onChunk(A2L_,rain)              = rain              (:)
    IF (PRESENT(snow))              dsl4jsb_var2D_onChunk(A2L_,snow)              = snow              (:)
    IF (PRESENT(wind_air))          dsl4jsb_var2D_onChunk(A2L_,wind_air)          = wind_air          (:)
    IF (PRESENT(wind_10m))          dsl4jsb_var2D_onChunk(A2L_,wind_10m)          = wind_10m          (:)
    IF (PRESENT(lw_srf_down))       dsl4jsb_var2D_onChunk(A2L_,lw_srf_down)       = lw_srf_down       (:)
    IF (PRESENT(swvis_srf_down))    dsl4jsb_var2D_onChunk(A2L_,swvis_srf_down)    = swvis_srf_down    (:)
    IF (PRESENT(swnir_srf_down))    dsl4jsb_var2D_onChunk(A2L_,swnir_srf_down)    = swnir_srf_down    (:)
    IF (PRESENT(swpar_srf_down))    dsl4jsb_var2D_onChunk(A2L_,swpar_srf_down)    = swpar_srf_down    (:)
! jsbach4    IF (PRESENT(frac_par_diffuse))  dsl4jsb_var2D_onChunk(A2L_,frac_par_diffuse)  = frac_par_diffuse  (:)
    IF (PRESENT(swpar_fdiffuse))    dsl4jsb_var2D_onChunk(A2L_,swpar_fdiffuse)    = swpar_fdiffuse    (:)
    IF (PRESENT(press_srf))         dsl4jsb_var2D_onChunk(A2L_,press_srf)         = press_srf         (:)
    IF (PRESENT(drag_srf))          dsl4jsb_var2D_onChunk(A2L_,drag_srf)          = drag_srf          (:)
    IF (PRESENT(t_acoef))           dsl4jsb_var2D_onChunk(A2L_,t_acoef)           = t_acoef           (:)
    IF (PRESENT(t_bcoef))           dsl4jsb_var2D_onChunk(A2L_,t_bcoef)           = t_bcoef           (:)
    IF (PRESENT(q_acoef))           dsl4jsb_var2D_onChunk(A2L_,q_acoef)           = q_acoef           (:)
    IF (PRESENT(q_bcoef))           dsl4jsb_var2D_onChunk(A2L_,q_bcoef)           = q_bcoef           (:)
    IF (PRESENT(pch))               dsl4jsb_var2D_onChunk(A2L_,pch)               = pch               (:)
    IF (PRESENT(cos_zenith_angle))  dsl4jsb_var2D_onChunk(A2L_,cos_zenith_angle)  = cos_zenith_angle  (:)
    IF (PRESENT(co2_mixing_ratio))     dsl4jsb_var2D_onChunk(A2L_,co2_mixing_ratio)     = co2_mixing_ratio     (:)
    IF (PRESENT(co2_mixing_ratio_C13)) dsl4jsb_var2D_onChunk(A2L_,co2_mixing_ratio_C13) = co2_mixing_ratio_C13 (:)
    IF (PRESENT(co2_mixing_ratio_C14)) dsl4jsb_var2D_onChunk(A2L_,co2_mixing_ratio_C14) = co2_mixing_ratio_C14 (:)
    IF (PRESENT(ga))                dsl4jsb_var2D_onChunk(A2L_,ga)                 = ga               (:)
    IF (PRESENT(nhx_deposition))    dsl4jsb_var2D_onChunk(A2L_,nhx_deposition)    = nhx_deposition    (:)
    IF (PRESENT(noy_deposition))    dsl4jsb_var2D_onChunk(A2L_,noy_deposition)    = noy_deposition    (:)
    IF (PRESENT(p_deposition))      dsl4jsb_var2D_onChunk(A2L_,p_deposition)      = p_deposition      (:)

    ! For lakes:
!    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      IF (PRESENT(drag_wtr))        dsl4jsb_var2D_onChunk(A2L_,drag_wtr)          = drag_wtr          (:)
      IF (PRESENT(drag_ice))        dsl4jsb_var2D_onChunk(A2L_,drag_ice)          = drag_ice          (:)
      IF (PRESENT(t_acoef_wtr))     dsl4jsb_var2D_onChunk(A2L_,t_acoef_wtr)       = t_acoef_wtr       (:)
      IF (PRESENT(t_bcoef_wtr))     dsl4jsb_var2D_onChunk(A2L_,t_bcoef_wtr)       = t_bcoef_wtr       (:)
      IF (PRESENT(q_acoef_wtr))     dsl4jsb_var2D_onChunk(A2L_,q_acoef_wtr)       = q_acoef_wtr       (:)
      IF (PRESENT(q_bcoef_wtr))     dsl4jsb_var2D_onChunk(A2L_,q_bcoef_wtr)       = q_bcoef_wtr       (:)
      IF (PRESENT(t_acoef_ice))     dsl4jsb_var2D_onChunk(A2L_,t_acoef_ice)       = t_acoef_ice       (:)
      IF (PRESENT(t_bcoef_ice))     dsl4jsb_var2D_onChunk(A2L_,t_bcoef_ice)       = t_bcoef_ice       (:)
      IF (PRESENT(q_acoef_ice))     dsl4jsb_var2D_onChunk(A2L_,q_acoef_ice)       = q_acoef_ice       (:)
      IF (PRESENT(q_bcoef_ice))     dsl4jsb_var2D_onChunk(A2L_,q_bcoef_ice)       = q_bcoef_ice       (:)
      IF (PRESENT(albedo_lwtr))     dsl4jsb_var2D_onChunk(A2L_,albedo_lwtr)       = albedo_lwtr       (:)
!    END IF

  END SUBROUTINE update_a2l


  !-----------------------------------------------------------------------------------------------------
  !> pass atmland variables from l2a memory back into the driver, per timestep
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_l2a( &
    ! the sequence of input variables is consistent with the sequence of var in t_l2a_memory !
    ! some variables differ compared to the JSBACH4 routine, see comments 
    ! all variables are defined as: real(wp), optional, intent(out)
      tile              ,   &
      options           ,   &
      t_srf             ,   &
      t_eff_srf         ,   &
      qsat_srf          ,   &
      fact_q_air        ,   &
      fact_qsat_srf     ,   &
      zh_srf            ,   &
      zm_srf            ,   &
      alb_vis_dif       ,   &
      alb_nir_dif       ,   &
      alb_vis_dir       ,   &
      alb_nir_dir       ,   &
      albedo            ,   &   ! added from SLM, not existent in jsbach4 20160823
      sensible_heat_flx ,   &   ! added from SLM, not existent in jsbach4 20160823
      latent_heat_flx   ,   &   ! added from SLM, not existent in jsbach4 20160823
      net_radiation     ,   &   ! added from SLM, not existent in jsbach4 20160823
      vpd               ,   &   ! added from SLM, not existent in jsbach4 20160823
      ! For lakes:
      latent_hflx_wtr   ,   &
      sensible_hflx_wtr ,   &
      latent_hflx_ice   ,   &
      sensible_hflx_ice     &
      & )

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< one tile
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    REAL(wp), OPTIONAL, INTENT(out) :: &
                                      t_srf(:)             ,   &
                                      t_eff_srf(:)         ,   &
                                      qsat_srf(:)          ,   &
                                      fact_q_air(:)        ,   &
                                      fact_qsat_srf(:)     ,   &
                                      zh_srf(:)            ,   &
                                      zm_srf(:)            ,   &
                                      alb_vis_dif(:)       ,   &
                                      alb_nir_dif(:)       ,   &
                                      alb_vis_dir(:)       ,   &
                                      alb_nir_dir(:)       ,   &
                                      albedo(:)            ,   &
                                      sensible_heat_flx(:) ,   &
                                      latent_heat_flx(:)   ,   &
                                      net_radiation(:)     ,   &
                                      vpd(:)               ,   &
                                      ! For lakes:
                                      latent_hflx_wtr(:)   ,   &
                                      sensible_hflx_wtr(:) ,   &
                                      latent_hflx_ice(:)   ,   &
                                      sensible_hflx_ice(:)
    ! ---------------------------
    ! 0.2 Local
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_l2a'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(L2A_)   !     TYPE(t_l2a_memory), POINTER :: l2a__mem
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(L2A_)) RETURN   ! ask for all processes considered here ?
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(L2A_)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! 1.0 pass variables values from memory to output variables
    IF (PRESENT(t_srf))                 t_srf(:)             = dsl4jsb_var2D_onChunk(L2A_,   t_srf)
    IF (PRESENT(t_eff_srf))             t_eff_srf(:)         = dsl4jsb_var2D_onChunk(L2A_,   t_eff_srf)
    IF (PRESENT(qsat_srf))              qsat_srf(:)          = dsl4jsb_var2D_onChunk(L2A_,   qsat_srf)
    IF (PRESENT(fact_q_air))            fact_q_air(:)        = dsl4jsb_var2D_onChunk(L2A_,   fact_q_air)
    IF (PRESENT(fact_qsat_srf))         fact_qsat_srf(:)     = dsl4jsb_var2D_onChunk(L2A_,   fact_qsat_srf)
    IF (PRESENT(zh_srf))                zh_srf(:)            = dsl4jsb_var2D_onChunk(L2A_,   zh_srf)
    IF (PRESENT(zm_srf))                zm_srf(:)            = dsl4jsb_var2D_onChunk(L2A_,   zm_srf)
    IF (PRESENT(alb_vis_dif))           alb_vis_dif(:)       = dsl4jsb_var2D_onChunk(L2A_,   alb_vis_dif)
    IF (PRESENT(alb_nir_dif))           alb_nir_dif(:)       = dsl4jsb_var2D_onChunk(L2A_,   alb_nir_dif)
    IF (PRESENT(alb_vis_dir))           alb_vis_dir(:)       = dsl4jsb_var2D_onChunk(L2A_,   alb_vis_dir)
    IF (PRESENT(alb_nir_dir))           alb_nir_dir(:)       = dsl4jsb_var2D_onChunk(L2A_,   alb_nir_dir)
    IF (PRESENT(albedo))                albedo(:)            = dsl4jsb_var2D_onChunk(L2A_,   albedo)
    IF (PRESENT(sensible_heat_flx))     sensible_heat_flx(:) = dsl4jsb_var2D_onChunk(L2A_,   sensible_heat_flx)
    IF (PRESENT(latent_heat_flx))       latent_heat_flx(:)   = dsl4jsb_var2D_onChunk(L2A_,   latent_heat_flx)
    IF (PRESENT(net_radiation))         net_radiation(:)     = dsl4jsb_var2D_onChunk(L2A_,   net_radiation)
    IF (PRESENT(vpd))                   vpd(:)               = dsl4jsb_var2D_onChunk(L2A_,   vpd)
    ! For lakes:
!    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      IF (PRESENT(latent_hflx_wtr))       latent_hflx_wtr(:)   = dsl4jsb_var2D_onChunk(L2A_,   latent_hflx_wtr)
      IF (PRESENT(sensible_hflx_wtr))     sensible_hflx_wtr(:) = dsl4jsb_var2D_onChunk(L2A_,   sensible_hflx_wtr)
      IF (PRESENT(latent_hflx_ice))       latent_hflx_ice(:)   = dsl4jsb_var2D_onChunk(L2A_,   latent_hflx_ice)
      IF (PRESENT(sensible_hflx_ice))     sensible_hflx_ice(:) = dsl4jsb_var2D_onChunk(L2A_,   sensible_hflx_ice)
!    ELSE
!      IF (PRESENT(latent_hflx_wtr))       latent_hflx_wtr(:)   = 0.0_wp
!      IF (PRESENT(sensible_hflx_wtr))     sensible_hflx_wtr(:) = 0.0_wp
!      IF (PRESENT(latent_hflx_ice))       latent_hflx_ice(:)   = 0.0_wp
!      IF (PRESENT(sensible_hflx_ice))     sensible_hflx_ice(:) = 0.0_wp
!    END IF

  END SUBROUTINE update_l2a


#endif
END MODULE mo_atmland_interface
