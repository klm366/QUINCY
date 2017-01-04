!>
!! @brief Contains the routines for the Task update_soil_physics
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Jan Engel (20170104): copyied from vegetation/mo_veg_soil.f90
!!
MODULE mo_soil_update_physics

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps8

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_soil_physics


  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_update_physics'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil moisture and theta for the plant only model
  !! 
  !! two cases considered: \n a) dynamic water uptake calculation, \n b) prescribed values from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_soil_physics(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,      ONLY: r_gas_dryair, LatentHeatVaporization, rhoh2o, grav, stbo, zemiss_def
    USE mo_phy_schemes
    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_veg_constants,         ONLY: k_soil_phi1, k_soil_phi2


    ! Use of process memories
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    dsl4jsb_Use_memory(LAND2ATM_)   !  USE mo_land2atm_memory_class, ONLY: t_land2atm_memory
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile     !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: t
    REAL(wp)                              :: q_sat,q_sat_soil
    REAL(wp)                              :: evaporation
    INTEGER                               :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_physics'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)                   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(LAND2ATM_)                   !     TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(SOIL_)                       !     TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_assimi_memory),   POINTER :: assimi          ! local pointer to assimi__mem for SLM usability
    TYPE(t_soil_memory),     POINTER :: soil             ! local pointer to veg__mem for SLM usability
    TYPE(t_veg_memory),      POINTER :: veg             ! local pointer to veg__mem for SLM usability
    REAL(wp), POINTER                :: t_air
    REAL(wp), POINTER                :: lw_srf_down
    REAL(wp), POINTER                :: press_srf
    REAL(wp), POINTER                :: q_air
    REAL(wp), POINTER                :: ga
    REAL(wp), POINTER                :: rain
    REAL(wp), POINTER                :: snow
    REAL(wp), POINTER                :: net_radiation
    REAL(wp), POINTER                :: ground_heat_flx
    REAL(wp), POINTER                :: latent_heat_flx
    REAL(wp), POINTER                :: sensible_heat_flx
    REAL(wp), POINTER                :: sw_srf_net
    REAL(wp), POINTER                :: drag_srf 
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(LAND2ATM_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(VEG_)

    assimi  => assimi__mem 
    soil    => soil__mem 
    veg     => veg__mem 
    t_air       => atm2land__mem%t_air%ptr
    lw_srf_down => atm2land__mem%lw_srf_down%ptr
    press_srf   => atm2land__mem%press_srf%ptr
    q_air       => atm2land__mem%q_air%ptr
    ga          => atm2land__mem%ga%ptr
    rain        => atm2land__mem%rain%ptr
    snow        => atm2land__mem%snow%ptr
    drag_srf    => atm2land__mem%drag_srf%ptr
    net_radiation      => land2atm__mem%net_radiation%ptr
    ground_heat_flx    => land2atm__mem%ground_heat_flx%ptr
    latent_heat_flx    => land2atm__mem%latent_heat_flx%ptr
    sensible_heat_flx  => land2atm__mem%sensible_heat_flx%ptr
    sw_srf_net         => rad__mem%sw_srf_net%ptr

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)


    ! saturating specific humidity 
    q_sat = calc_spec_humidity_sat(t_air,press_srf)
    q_sat_soil = calc_spec_humidity_sat(soil%t_soil%ptr(1,1,1),press_srf)

    IF(assimi%canopy_cond%ptr > eps8)THEN
       ! T = rho * (qsat - qair ) * Gtot in mm/s 
       veg%transpiration%ptr = press_srf / ( r_gas_dryair * t_air ) * & 
                  (q_sat - q_air) / (1._wp / assimi%canopy_cond%ptr + 1._wp / ga )  
    ELSE 
       veg%transpiration%ptr = 0.0_wp 
    ENDIF
    ! E = rho * (qsat - qair) * theta * 1/ga_soil (taken as ga_soil = ga / 10)
    evaporation = press_srf / ( r_gas_dryair * t_air ) * & 
                  (q_sat - q_air) / ( 100._wp / ga) * & 
                  soil%w_soil_root_theta%ptr(1,1) 
    IF(sw_srf_net < eps8) evaporation = 0.0_wp

    ! calculate latent heat flux (W m-2) 
    latent_heat_flx = LatentHeatVaporization * (evaporation + veg%transpiration%ptr)

    ! estimate surface temperature and sensible heat flux 
    IF (soil_config%flag_calc_soil_energy) THEN ! case of a simple diagnostic soil energy calculation
       ! simple estimate of sensible heat flux (ignoring instabilities and actual surface moisture...)
       sensible_heat_flx = heat_transfer_coef(drag_srf, model%options%steplen) * & 
                             (surface_dry_static_energy(soil%t_soil%ptr(1,1,1),q_sat_soil) - & 
                              surface_dry_static_energy(t_air,q_air))
       ! update deep soil temperature to estimate ground flux from the surface layer
       ! ignores interaction of diffusion with water or ice
       DO isoil = soil_config%nsoil_water,2,-1
          IF(soil%soil_depth_l%ptr(1,isoil,1)>eps8)THEN
             ground_heat_flx = (soil%t_soil%ptr(1,isoil-1,1)-soil%t_soil%ptr(1,isoil,1)) * soil%therm_cond_l%ptr(1,isoil,1) / & 
                        (soil%soil_depth_l%ptr(1,isoil,1)/2.0_wp+soil%soil_depth_l%ptr(1,isoil-1,1)/2.0_wp)
             soil%t_soil%ptr(:,isoil,1) = soil%t_soil%ptr(:,isoil,1) + ground_heat_flx * dtime / & 
                        (soil%heat_capa_l%ptr(:,isoil,1) * soil%bulk_dens_l%ptr(:,isoil,1) * soil%soil_depth_l%ptr(:,isoil,1)) 
          ENDIF 
       ENDDO
       ! update net radiation flux at the surface
       net_radiation = latent_heat_flx + sensible_heat_flx - ground_heat_flx 
       soil%t_soil%ptr(:,1,1) = soil%t_soil%ptr(:,1,1) + ( sw_srf_net - net_radiation + & 
                                ! simple estimate of surface net longwave flux
                                lw_srf_down - stbo * zemiss_def * soil%t_soil%ptr(:,1,1)**4 ) / & 
                                (soil%heat_capa_l%ptr(:,1,1) * soil%bulk_dens_l%ptr(:,1,1) * soil%soil_depth_l%ptr(:,1,1)) * &
                                dtime 

    ELSE ! diagnose sensible heat flux from net radiation and set soil temperature to air temperature 
       ground_heat_flx = 0.0_wp 
       sensible_heat_flx = net_radiation - latent_heat_flx 
       soil%t_soil%ptr(:,:,1) = t_air
    ENDIF

   IF(soil_config%flag_calc_soil_moisture) THEN ! case of dynamic water uptake calculation 

       ! update soil moisture dS = rain - transpiration - evaporation, incl. conversion to [m]
       soil%w_soil_root%ptr(1,1) = soil%w_soil_root%ptr(1,1) + (rain + snow &
            - veg%transpiration%ptr - evaporation) / 1000._wp * dtime

       ! calculate drainage if any [mm/s]
       WHERE(soil%w_soil_root%ptr > soil%w_soil_root_fc%ptr)
         soil%drainage%ptr = ( soil%w_soil_root%ptr - soil%w_soil_root_fc%ptr ) * 1000._wp/ dtime
       ELSEWHERE
         soil%drainage%ptr = 0.0_wp
       ENDWHERE
       soil%w_soil_root%ptr = soil%w_soil_root%ptr - soil%drainage%ptr / 1000.0_wp * dtime
       soil%w_soil_root%ptr = MAX(soil%w_soil_root%ptr,0.0_wp)

       ! scale to relative water availability (fraction of maximum)
       soil%w_soil_root_theta%ptr = ( soil%w_soil_root%ptr - soil%w_soil_root_pwp%ptr ) / &
                                    ( soil%w_soil_root_fc%ptr - soil%w_soil_root_pwp%ptr ) 
       soil%w_soil_root_theta%ptr = MAX(soil%w_soil_root_theta%ptr,0.0_wp)
  
   ENDIF 

   ! calculate soil wtaer potential according to Campbell et al 1974 
   soil%w_soil_root_pot%ptr = k_soil_phi1 * & 
         ( soil%w_soil_root%ptr / soil%w_soil_root_fc%ptr )**k_soil_phi2

   veg%dphi%ptr = soil%w_soil_root_pot%ptr(1,1) - lctlib%phi_leaf_min - grav * rhoh2o * veg%height%ptr * 1.e-6_wp

  END SUBROUTINE update_soil_physics



END MODULE mo_soil_update_physics
