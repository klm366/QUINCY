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
  USE mo_base_memory,             ONLY: t_slm_tile
  USE mo_jsb_var_class,           ONLY: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_real1d
  USE mo_jsb_process_class,       ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_

  ! Use of process memories
  dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(LAND2ATM_)   !  USE mo_land2atm_memory_class, ONLY: t_land2atm_memory

!!! jsbach4
!!  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
!!  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
!!  USE mo_jsb_config_class,   ONLY: t_jsb_config
!!  USE mo_jsb_process_class,  ONLY: t_jsb_process
!!  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
!!  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE, Contains_lct

!!  dsl4jsb_Use_processes ATM2LAND_, SEB_, TURB_, SOIL_, HYDRO_, RAD_
!!  dsl4jsb_Use_memory(ATM2LAND_)
!!  dsl4jsb_Use_memory(SEB_)
!!  dsl4jsb_Use_memory(TURB_)
!!  dsl4jsb_Use_memory(SOIL_)
!!  dsl4jsb_Use_memory(HYDRO_)
!!  dsl4jsb_Use_memory(RAD_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_atm2land, update_land2atm


  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> pass atmland variables from driver into atm2land memory per timestep
  !!
  !! Routine is preliminary, as it is adapted to work with current SLM configurations (e.g. REAL but 2D)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_atm2land( &
    ! the sequence of input variables is consistent with the sequence of var in t_atm2land_memory !
    ! some variables differ compared to the JSBACH4 routine, see comments 
    ! all variables are defined as: real(wp), optional, intent(in)
      tile               , &
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
      ga                   , &  !> @todo ga in t_atm2land_memory: is a variable from SLM, which may be moved elsewhere
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


    IMPLICIT NONE 

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout) :: tile           !< one tile
    REAL(wp), OPTIONAL, INTENT(in) :: &
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
                                      ga                   , &  !> @todo ga in t_atm2land_memory: is a variable from SLM, which may be moved elsewhere
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
                                      q_bcoef_ice
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_atm2land'

    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Get_memory(ATM2LAND_)
!    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p) 
!      TYPE IS(t_atm2land_memory)
!        atm2land__mem => mm
!    END SELECT


    ! 1.0 pass variables values to memory
    !IF (PRESENT(t_air))             dsl4jsb_var2D_onChunk(ATM2LAND_,t_air)             = t_air             (:)
    !from dsl_jsbach.inc: #define dsl4jsb_var2D_onChunk(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,iblk)
    IF (PRESENT(t_air))                   atm2land__mem%t_air%ptr                       = t_air
    IF (PRESENT(q_air))                   atm2land__mem%q_air%ptr                       = q_air
    IF (PRESENT(rain))                    atm2land__mem%rain%ptr                        = rain
    IF (PRESENT(snow))                    atm2land__mem%snow%ptr                        = snow
    IF (PRESENT(wind_air))                atm2land__mem%wind_air%ptr                    = wind_air
    IF (PRESENT(wind_10m))                atm2land__mem%wind_10m%ptr                    = wind_10m
    IF (PRESENT(lw_srf_down))             atm2land__mem%lw_srf_down%ptr                 = lw_srf_down
    IF (PRESENT(swvis_srf_down))          atm2land__mem%swvis_srf_down%ptr              = swvis_srf_down     
    IF (PRESENT(swnir_srf_down))          atm2land__mem%swnir_srf_down%ptr              = swnir_srf_down     
    IF (PRESENT(swpar_srf_down))          atm2land__mem%swpar_srf_down%ptr              = swpar_srf_down     
    IF (PRESENT(swpar_fdiffuse))          atm2land__mem%swpar_fdiffuse%ptr              = swpar_fdiffuse     
    IF (PRESENT(press_srf))               atm2land__mem%press_srf%ptr                   = press_srf
    IF (PRESENT(drag_srf))                atm2land__mem%drag_srf%ptr                    = drag_srf
    IF (PRESENT(t_acoef))                 atm2land__mem%t_acoef%ptr                     = t_acoef
    IF (PRESENT(t_bcoef))                 atm2land__mem%t_bcoef%ptr                     = t_bcoef
    IF (PRESENT(q_acoef))                 atm2land__mem%q_acoef%ptr                     = q_acoef
    IF (PRESENT(q_bcoef))                 atm2land__mem%q_bcoef%ptr                     = q_bcoef
    IF (PRESENT(pch))                     atm2land__mem%pch%ptr                         = pch
    IF (PRESENT(cos_zenith_angle))        atm2land__mem%cos_zenith_angle%ptr            = cos_zenith_angle   
    IF (PRESENT(co2_mixing_ratio))        atm2land__mem%co2_mixing_ratio%ptr            = co2_mixing_ratio     
    IF (PRESENT(co2_mixing_ratio_C13))    atm2land__mem%co2_mixing_ratio_C13%ptr        = co2_mixing_ratio_C13
    IF (PRESENT(co2_mixing_ratio_C14))    atm2land__mem%co2_mixing_ratio_C14%ptr        = co2_mixing_ratio_C14
    IF (PRESENT(ga))                      atm2land__mem%ga%ptr                          = ga
    ! For lakes:
    IF (PRESENT(albedo_lwtr))             atm2land__mem%albedo_lwtr%ptr                 = albedo_lwtr
    IF (PRESENT(drag_wtr))                atm2land__mem%drag_wtr%ptr                    = drag_wtr
    IF (PRESENT(drag_ice))                atm2land__mem%drag_ice%ptr                    = drag_ice
    IF (PRESENT(t_acoef_wtr))             atm2land__mem%t_acoef_wtr%ptr                 = t_acoef_wtr
    IF (PRESENT(t_bcoef_wtr))             atm2land__mem%t_bcoef_wtr%ptr                 = t_bcoef_wtr
    IF (PRESENT(q_acoef_wtr))             atm2land__mem%q_acoef_wtr%ptr                 = q_acoef_wtr
    IF (PRESENT(q_bcoef_wtr))             atm2land__mem%q_bcoef_wtr%ptr                 = q_bcoef_wtr
    IF (PRESENT(t_acoef_ice))             atm2land__mem%t_acoef_ice%ptr                 = t_acoef_ice
    IF (PRESENT(t_bcoef_ice))             atm2land__mem%t_bcoef_ice%ptr                 = t_bcoef_ice
    IF (PRESENT(q_acoef_ice))             atm2land__mem%q_acoef_ice%ptr                 = q_acoef_ice
    IF (PRESENT(q_bcoef_ice))             atm2land__mem%q_bcoef_ice%ptr                 = q_bcoef_ice

  END SUBROUTINE update_atm2land
!!! jsbach4
!!  SUBROUTINE update_atm2land( &
!!    & tile, options,          &
!!    & t_air,                  &
!!    & q_air,                  &
!!    & rain,                   &
!!    & snow,                   &
!!    & wind_air,               & 
!!    & wind_10m,               &
!!    & lw_srf_down,            &
!!    & swvis_srf_down,         &
!!    & swnir_srf_down,         &
!!    & swpar_srf_down,         &
!!    & frac_par_diffuse,       &
!!    & press_srf,              &
!!    & drag_srf,               &
!!    & t_acoef,                &
!!    & t_bcoef,                &
!!    & q_acoef,                &
!!    & q_bcoef,                &
!!    & pch,                    &
!!    & cos_zenith_angle,       &
!!    & CO2_air,                &
!!    & drag_wtr,               &
!!    & drag_ice,               &
!!    & t_acoef_wtr,            &
!!    & t_bcoef_wtr,            &
!!    & q_acoef_wtr,            &
!!    & q_bcoef_wtr,            &
!!    & t_acoef_ice,            &
!!    & t_bcoef_ice,            &
!!    & q_acoef_ice,            &
!!    & q_bcoef_ice,            &
!!    & albedo_lwtr             &
!!    & )

!!    USE mo_jsb_physical_constants, ONLY: molarMassDryAir, molarMassCO2

!!    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
!!    TYPE(t_jsb_task_options),   INTENT(in)    :: options
!!    REAL(wp), OPTIONAL, INTENT(in) ::                                            &
!!      & t_air(:), q_air(:), rain(:), snow(:), wind_air(:), wind_10m(:),          &
!!      & lw_srf_down(:), swvis_srf_down(:), swnir_srf_down(:), swpar_srf_down(:), &
!!      & frac_par_diffuse(:), press_srf(:), drag_srf(:),                                               &
!!      & t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:),                          &
!!      & pch(:), cos_zenith_angle(:),                                             &
!!      & CO2_air(:),                                                              &
!!      & drag_wtr(:), drag_ice(:),                                                &
!!      & t_acoef_wtr(:), t_bcoef_wtr(:), q_acoef_wtr(:), q_bcoef_wtr(:),          &
!!      & t_acoef_ice(:), t_bcoef_ice(:), q_acoef_ice(:), q_bcoef_ice(:),          &
!!      & albedo_lwtr(:)

!!    dsl4jsb_Def_memory(ATM2LAND_)

!!    INTEGER ::  iblk, ics, ice

!!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_atm2land'

!!    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

!!    iblk = options%iblk
!!    ics  = options%ics
!!    ice  = options%ice

!!    IF (ice-ics+1 /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

!!    dsl4jsb_Get_memory(ATM2LAND_)

!!    IF (PRESENT(t_air))             dsl4jsb_var2D_onChunk(ATM2LAND_,t_air)             = t_air             (:)
!!    IF (PRESENT(q_air))             dsl4jsb_var2D_onChunk(ATM2LAND_,q_air)             = q_air             (:)
!!    IF (PRESENT(rain))              dsl4jsb_var2D_onChunk(ATM2LAND_,rain)              = rain              (:)
!!    IF (PRESENT(snow))              dsl4jsb_var2D_onChunk(ATM2LAND_,snow)              = snow              (:)
!!    IF (PRESENT(wind_air))          dsl4jsb_var2D_onChunk(ATM2LAND_,wind_air)          = wind_air          (:)
!!    IF (PRESENT(wind_10m))          dsl4jsb_var2D_onChunk(ATM2LAND_,wind_10m)          = wind_10m          (:)
!!    IF (PRESENT(lw_srf_down))       dsl4jsb_var2D_onChunk(ATM2LAND_,lw_srf_down)       = lw_srf_down       (:)
!!    IF (PRESENT(swvis_srf_down))    dsl4jsb_var2D_onChunk(ATM2LAND_,swvis_srf_down)    = swvis_srf_down    (:)
!!    IF (PRESENT(swnir_srf_down))    dsl4jsb_var2D_onChunk(ATM2LAND_,swnir_srf_down)    = swnir_srf_down    (:)
!!    IF (PRESENT(swpar_srf_down))    dsl4jsb_var2D_onChunk(ATM2LAND_,swpar_srf_down)    = swpar_srf_down    (:)
!!    IF (PRESENT(frac_par_diffuse))  dsl4jsb_var2D_onChunk(ATM2LAND_,frac_par_diffuse)  = frac_par_diffuse  (:)
!!    IF (PRESENT(press_srf))         dsl4jsb_var2D_onChunk(ATM2LAND_,press_srf)         = press_srf         (:)
!!    IF (PRESENT(drag_srf))          dsl4jsb_var2D_onChunk(ATM2LAND_,drag_srf)          = drag_srf          (:)
!!    IF (PRESENT(t_acoef))           dsl4jsb_var2D_onChunk(ATM2LAND_,t_acoef)           = t_acoef           (:)
!!    IF (PRESENT(t_bcoef))           dsl4jsb_var2D_onChunk(ATM2LAND_,t_bcoef)           = t_bcoef           (:)
!!    IF (PRESENT(q_acoef))           dsl4jsb_var2D_onChunk(ATM2LAND_,q_acoef)           = q_acoef           (:)
!!    IF (PRESENT(q_bcoef))           dsl4jsb_var2D_onChunk(ATM2LAND_,q_bcoef)           = q_bcoef           (:)
!!    IF (PRESENT(pch))               dsl4jsb_var2D_onChunk(ATM2LAND_,pch)               = pch               (:)
!!    IF (PRESENT(cos_zenith_angle))  dsl4jsb_var2D_onChunk(ATM2LAND_,cos_zenith_angle)  = cos_zenith_angle  (:)
!!    IF (PRESENT(CO2_air)) THEN
!!      dsl4jsb_var2D_onChunk(ATM2LAND_, CO2_air)               = CO2_air(:)
!!      ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
!!      dsl4jsb_var2D_onChunk(ATM2LAND_, CO2_air_mol) = CO2_air(:) * molarMassDryAir / molarMassCO2
!!    END IF

!!    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN

!!      IF (PRESENT(drag_wtr))        dsl4jsb_var2D_onChunk(ATM2LAND_,drag_wtr)          = drag_wtr          (:)
!!      IF (PRESENT(drag_ice))        dsl4jsb_var2D_onChunk(ATM2LAND_,drag_ice)          = drag_ice          (:)
!!      IF (PRESENT(t_acoef_wtr))     dsl4jsb_var2D_onChunk(ATM2LAND_,t_acoef_wtr)       = t_acoef_wtr       (:)
!!      IF (PRESENT(t_bcoef_wtr))     dsl4jsb_var2D_onChunk(ATM2LAND_,t_bcoef_wtr)       = t_bcoef_wtr       (:)
!!      IF (PRESENT(q_acoef_wtr))     dsl4jsb_var2D_onChunk(ATM2LAND_,q_acoef_wtr)       = q_acoef_wtr       (:)
!!      IF (PRESENT(q_bcoef_wtr))     dsl4jsb_var2D_onChunk(ATM2LAND_,q_bcoef_wtr)       = q_bcoef_wtr       (:)
!!      IF (PRESENT(t_acoef_ice))     dsl4jsb_var2D_onChunk(ATM2LAND_,t_acoef_ice)       = t_acoef_ice       (:)
!!      IF (PRESENT(t_bcoef_ice))     dsl4jsb_var2D_onChunk(ATM2LAND_,t_bcoef_ice)       = t_bcoef_ice       (:)
!!      IF (PRESENT(q_acoef_ice))     dsl4jsb_var2D_onChunk(ATM2LAND_,q_acoef_ice)       = q_acoef_ice       (:)
!!      IF (PRESENT(q_bcoef_ice))     dsl4jsb_var2D_onChunk(ATM2LAND_,q_bcoef_ice)       = q_bcoef_ice       (:)
!!      IF (PRESENT(albedo_lwtr))     dsl4jsb_var2D_onChunk(ATM2LAND_,albedo_lwtr)       = albedo_lwtr       (:)

!!    END IF

!!  END SUBROUTINE update_atm2land


  !-----------------------------------------------------------------------------------------------------
  !> pass atmland variables from land2atm memory back into the driver, per timestep
  !!
  !! Routine is preliminary, as it is adapted to work with current SLM configurations (e.g. REAL but 2D)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_land2atm( &
    ! the sequence of input variables is consistent with the sequence of var in t_land2atm_memory !
    ! some variables differ compared to the JSBACH4 routine, see comments 
    ! all variables are defined as: real(wp), optional, intent(out)
      tile              ,   &
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

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout) :: tile           !< one tile
    REAL(wp), OPTIONAL, INTENT(out) :: &
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
                                      albedo            ,   &
                                      sensible_heat_flx ,   &
                                      latent_heat_flx   ,   &
                                      net_radiation     ,   &
                                      vpd               ,   &
                                      ! For lakes:
                                      latent_hflx_wtr   ,   &
                                      sensible_hflx_wtr ,   &
                                      latent_hflx_ice   ,   &
                                      sensible_hflx_ice
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_land2atm'

    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(LAND2ATM_)   !     TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    dsl4jsb_Get_memory(LAND2ATM_)
!    SELECT TYPE (mm => tile%mem(LAND2ATM_)%p) 
!      TYPE IS(t_land2atm_memory)
!        land2atm__mem => mm
!    END SELECT


    ! 1.0 pass variables values from memory to output variables
    !IF (PRESENT(t_srf))           t_srf(:)            = dsl4jsb_var2D_onChunk(SEB_,   t_filt)
    !from dsl_jsbach.inc: #define dsl4jsb_var2D_onChunk(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,iblk)
    IF (PRESENT(t_srf))                 t_srf             = land2atm__mem%t_srf%ptr
    IF (PRESENT(t_eff_srf))             t_eff_srf         = land2atm__mem%t_eff_srf%ptr
    IF (PRESENT(qsat_srf))              qsat_srf          = land2atm__mem%qsat_srf%ptr
    IF (PRESENT(fact_q_air))            fact_q_air        = land2atm__mem%fact_q_air%ptr
    IF (PRESENT(fact_qsat_srf))         fact_qsat_srf     = land2atm__mem%fact_qsat_srf%ptr
    IF (PRESENT(zh_srf))                zh_srf            = land2atm__mem%zh_srf%ptr
    IF (PRESENT(zm_srf))                zm_srf            = land2atm__mem%zm_srf%ptr
    IF (PRESENT(alb_vis_dif))           alb_vis_dif       = land2atm__mem%alb_vis_dif%ptr
    IF (PRESENT(alb_nir_dif))           alb_nir_dif       = land2atm__mem%alb_nir_dif%ptr
    IF (PRESENT(alb_vis_dir))           alb_vis_dir       = land2atm__mem%alb_vis_dir%ptr
    IF (PRESENT(alb_nir_dir))           alb_nir_dir       = land2atm__mem%alb_nir_dir%ptr
    IF (PRESENT(albedo))                albedo            = land2atm__mem%albedo%ptr
    IF (PRESENT(sensible_heat_flx))     sensible_heat_flx = land2atm__mem%sensible_heat_flx%ptr
    IF (PRESENT(latent_heat_flx))       latent_heat_flx   = land2atm__mem%latent_heat_flx%ptr
    IF (PRESENT(net_radiation))         net_radiation     = land2atm__mem%net_radiation%ptr
    IF (PRESENT(vpd))                   vpd               = land2atm__mem%vpd%ptr
    ! For lakes:
    IF (PRESENT(latent_hflx_wtr))       latent_hflx_wtr   = land2atm__mem%latent_hflx_wtr%ptr
    IF (PRESENT(sensible_hflx_wtr))     sensible_hflx_wtr = land2atm__mem%sensible_hflx_wtr%ptr
    IF (PRESENT(latent_hflx_ice))       latent_hflx_ice   = land2atm__mem%latent_hflx_ice%ptr
    IF (PRESENT(sensible_hflx_ice))     sensible_hflx_ice = land2atm__mem%sensible_hflx_ice%ptr

  END SUBROUTINE update_land2atm
!!! jsbach4
!!  SUBROUTINE update_land2atm(tile, options,                                          &
!!    & t_srf, t_eff_srf, qsat_srf, s_srf,                                             &
!!    & fact_q_air, fact_qsat_srf, evapopot,                                           &
!!    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,                  &
!!    & rough_h_srf, rough_m_srf, tte_corr,                                            &
!!    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,                            &
!!    & t_lwtr, t_lice,                                                                &
!!    & evapo_wtr, latent_hflx_wtr, sensible_hflx_wtr,                                 &
!!    & evapo_ice, latent_hflx_ice, sensible_hflx_ice,                                 &
!!    & ice_fract_lake,                                                                &
!!    & alb_vis_dir_wtr, alb_vis_dif_wtr, alb_nir_dir_wtr, alb_nir_dif_wtr, albedo_lice)

!!    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
!!    TYPE(t_jsb_task_options),   INTENT(in)    :: options
!!    REAL(wp), OPTIONAL, INTENT(out) ::                                               &
!!      & t_srf(:), t_eff_srf(:), qsat_srf(:), s_srf(:),                               &
!!      & fact_q_air(:), fact_qsat_srf(:), evapopot(:),                                &
!!      & evapotrans(:), latent_hflx(:), sensible_hflx(:), grnd_hflx(:), grnd_hcap(:), &
!!      & rough_h_srf(:), rough_m_srf(:), tte_corr(:),                                 &
!!      & t_lwtr(:), t_lice(:),                                                        &
!!      & evapo_wtr(:), latent_hflx_wtr(:), sensible_hflx_wtr(:),                      &
!!      & evapo_ice(:), latent_hflx_ice(:), sensible_hflx_ice(:),                      &
!!      & ice_fract_lake(:),                                                           &
!!      & alb_vis_dir_wtr(:), alb_vis_dif_wtr(:), alb_nir_dir_wtr(:),                  &
!!      & alb_nir_dif_wtr(:), albedo_lice(:),                                             &
!!      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:)

!!    dsl4jsb_Def_memory(SEB_)
!!    dsl4jsb_Def_memory(TURB_)
!!    dsl4jsb_Def_memory(SOIL_)
!!    dsl4jsb_Def_memory(HYDRO_)
!!    dsl4jsb_Def_memory(RAD_)

!!    INTEGER ::  iblk, ics, ice

!!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_land2atm'

!!    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

!!    iblk = options%iblk
!!    ics  = options%ics
!!    ice  = options%ice

!!    IF (ice-ics+1 /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

!!    dsl4jsb_Get_memory(SEB_)
!!    dsl4jsb_Get_memory(TURB_)
!!    dsl4jsb_Get_memory(SOIL_)
!!    dsl4jsb_Get_memory(HYDRO_)
!!    dsl4jsb_Get_memory(RAD_)

!!    IF (PRESENT(t_srf))           t_srf(:)            = dsl4jsb_var2D_onChunk(SEB_,   t_filt)
!!    IF (PRESENT(t_eff_srf))       t_eff_srf(:)        = dsl4jsb_var2D_onChunk(SEB_,   t_eff)
!!    IF (PRESENT(qsat_srf))        qsat_srf(:)         = dsl4jsb_var2D_onChunk(SEB_,   qsat_star)
!!    IF (PRESENT(s_srf))           s_srf(:)            = dsl4jsb_var2D_onChunk(SEB_,   s_star)
!!    IF (PRESENT(fact_q_air))      fact_q_air(:)       = dsl4jsb_var2D_onChunk(TURB_,  fact_q_air)
!!    IF (PRESENT(fact_qsat_srf))   fact_qsat_srf(:)    = dsl4jsb_var2D_onChunk(TURB_,  fact_qsat_srf)
!!    IF (PRESENT(evapopot))        evapopot(:)         = dsl4jsb_var2D_onChunk(HYDRO_, evapopot)
!!    IF (PRESENT(evapotrans))      evapotrans(:)       = dsl4jsb_var2D_onChunk(HYDRO_, evapo_lnd)
!!    IF (PRESENT(latent_hflx))     latent_hflx(:)      = dsl4jsb_var2D_onChunk(SEB_,   latent_hflx)
!!    IF (PRESENT(sensible_hflx))   sensible_hflx(:)    = dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx)
!!    IF (PRESENT(grnd_hflx))       grnd_hflx(:)        = dsl4jsb_var2D_onChunk(SEB_,   forc_hflx)
!!    IF (PRESENT(grnd_hcap))       grnd_hcap(:)        = dsl4jsb_var2D_onChunk(SEB_,   heat_cap_old)
!!    IF (PRESENT(rough_h_srf))     rough_h_srf(:)      = dsl4jsb_var2D_onChunk(TURB_,  rough_h)
!!    IF (PRESENT(rough_m_srf))     rough_m_srf(:)      = dsl4jsb_var2D_onChunk(TURB_,  rough_m)
!!    IF (PRESENT(tte_corr))        tte_corr(:)         = dsl4jsb_var2D_onChunk(HYDRO_, tte_corr)
!!    IF (PRESENT(alb_vis_dir))     alb_vis_dir(:)      = dsl4jsb_var2D_onChunk(RAD_,   alb_vis)  ! Temporary
!!    IF (PRESENT(alb_vis_dif))     alb_vis_dif(:)      = dsl4jsb_var2D_onChunk(RAD_,   alb_vis)  ! Temporary
!!    IF (PRESENT(alb_nir_dir))     alb_nir_dir(:)      = dsl4jsb_var2D_onChunk(RAD_,   alb_nir)  ! Temporary
!!    IF (PRESENT(alb_nir_dif))     alb_nir_dif(:)      = dsl4jsb_var2D_onChunk(RAD_,   alb_nir)  ! Temporary

!!    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
!!      IF (PRESENT(t_lwtr))            t_lwtr(:)            = dsl4jsb_var2D_onChunk(SEB_,   t_lwtr)
!!      IF (PRESENT(t_lice))            t_lice(:)            = dsl4jsb_var2D_onChunk(SEB_,   t_lice)
!!      IF (PRESENT(evapo_wtr))         evapo_wtr(:)         = dsl4jsb_var2D_onChunk(HYDRO_, evapo_wtr)
!!      IF (PRESENT(latent_hflx_wtr))   latent_hflx_wtr(:)   = dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_wtr)
!!      IF (PRESENT(sensible_hflx_wtr)) sensible_hflx_wtr(:) = dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_wtr)
!!      IF (PRESENT(evapo_ice))         evapo_ice(:)         = dsl4jsb_var2D_onChunk(HYDRO_, evapo_ice)
!!      IF (PRESENT(latent_hflx_ice))   latent_hflx_ice(:)   = dsl4jsb_var2D_onChunk(SEB_,   latent_hflx_ice)
!!      IF (PRESENT(sensible_hflx_ice)) sensible_hflx_ice(:) = dsl4jsb_var2D_onChunk(SEB_,   sensible_hflx_ice)
!!      IF (PRESENT(ice_fract_lake))    ice_fract_lake(:)    = dsl4jsb_var2D_onChunk(SEB_,   fract_lice)
!!      IF (PRESENT(albedo_lice))       albedo_lice(:)       = dsl4jsb_var2D_onChunk(RAD_,   albedo_lice)
!!    ELSE
!!      IF (PRESENT(t_lwtr))            t_lwtr(:)             = 0._wp                                   ! Temporary
!!      IF (PRESENT(t_lice))            t_lice(:)             = 0._wp                                   ! Temporary
!!      IF (PRESENT(evapo_wtr))         evapo_wtr(:)         = 0._wp
!!      IF (PRESENT(latent_hflx_wtr))   latent_hflx_wtr(:)   = 0._wp
!!      IF (PRESENT(sensible_hflx_wtr)) sensible_hflx_wtr(:) = 0._wp
!!      IF (PRESENT(evapo_ice))         evapo_ice(:)         = 0._wp
!!      IF (PRESENT(latent_hflx_ice))   latent_hflx_ice(:)   = 0._wp
!!      IF (PRESENT(sensible_hflx_ice)) sensible_hflx_ice(:) = 0._wp
!!      IF (PRESENT(ice_fract_lake))    ice_fract_lake(:)    = 0._wp                                   ! Temporary
!!      IF (PRESENT(albedo_lice))       albedo_lice(:)       = 0.55_wp
!!    END IF

!!  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
