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
    USE mo_soil_constants,        ONLY: k_soil_phi1, k_soil_phi2, drain_frac, fdrain_runoff, kdiff_theta, min_theta_wpot
    USE mo_assimi_update_canopy_fluxes, ONLY: calc_soil_moisture_stress

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
    REAL(wp)                              :: q_sat,q_sat_soil
    REAL(wp)                              :: evaporation
    REAL(wp)                              :: hlp1,infiltration
    REAL(wp),DIMENSION(soil_config%nsoil_water) :: ftranspiration,fevaporation,theta_sl
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
    TYPE(t_soil_memory),     POINTER :: soil            ! local pointer to soil__mem for SLM usability
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


    !> 1.0 Latent heat flux
    !!

    !> 1.1 calculate saturating specific humidity 
    !!
    q_sat = calc_spec_humidity_sat(t_air,press_srf)
    q_sat_soil = calc_spec_humidity_sat(soil%t_soil%ptr(1,1,1),press_srf)

    !> 1.2. calulate transpiration (kg / m2 / s)
    !! T = rho * (qsat - qair ) * Gtot in mm/s 
    !!
    IF(assimi%canopy_cond%ptr > eps8)THEN
       veg%transpiration%ptr = press_srf / ( r_gas_dryair * t_air ) * & 
                  (q_sat - q_air) / (1._wp / assimi%canopy_cond%ptr + 1._wp / ga )  
    ELSE 
       veg%transpiration%ptr = 0.0_wp 
    ENDIF

    !> 1.3 soil evaporation (kg / m2 / s)
    !! E = rho * (qsat - qair) * theta * 1/ga_soil (taken as ga_soil = ga / 10)
    !!
    SELECT CASE(soil_config%soil_moisture_scheme)
       CASE("bucket") ! ! Case of bucket water scheme
          evaporation = press_srf / ( r_gas_dryair * t_air ) * & 
                        (q_sat - q_air) / ( 100._wp / ga) * & 
                        soil%w_soil_root_theta%ptr(1,1) 
          IF(sw_srf_net < eps8) evaporation = 0.0_wp
       
       CASE("layer") ! Case of layer scheme
          evaporation = press_srf / ( r_gas_dryair * t_air ) * & 
                        (q_sat - q_air) / ( 10._wp / ga) * & 
                        SUM(soil%w_soil_sl%ptr(1,1:2,1))/SUM(soil%w_soil_fc_sl%ptr(1,1:2,1)) 
          IF(evaporation < eps8) evaporation = 0.0_wp
    ENDSELECT

    !> 1.4 calculate latent heat flux (W m-2) 
    !!
    latent_heat_flx = LatentHeatVaporization * (evaporation + veg%transpiration%ptr)

    !> 2.0 Estimate surface temperature and sensible heat flux 
    !!
    IF (soil_config%flag_calc_soil_energy) THEN ! case of a simple diagnostic soil energy calculation
       !> 2.1 simple estimate of sensible heat flux
       !! (ignoring instabilities and actual surface moisture...)
       !!
       sensible_heat_flx = heat_transfer_coef(drag_srf, model%options%steplen) * & 
                             (surface_dry_static_energy(soil%t_soil%ptr(1,1,1),q_sat_soil) - & 
                              surface_dry_static_energy(t_air,q_air))

       !> 2.2 Estimate ground flux and update deep soil temperature
       !! working upwards to estimate net transfer implied by previous temperature profile 
       !! (ignores interaction of diffusion with water or ice)
       !!
       DO isoil = soil_config%nsoil_water,2,-1
          IF(soil%soil_depth_l%ptr(1,isoil,1)>eps8)THEN
             ground_heat_flx = (soil%t_soil%ptr(1,isoil-1,1)-soil%t_soil%ptr(1,isoil,1)) * soil%therm_cond_l%ptr(1,isoil,1) / & 
                        (soil%soil_depth_l%ptr(1,isoil,1)/2.0_wp+soil%soil_depth_l%ptr(1,isoil-1,1)/2.0_wp)
             soil%t_soil%ptr(:,isoil,1) = soil%t_soil%ptr(:,isoil,1) + ground_heat_flx * dtime / & 
                        (soil%heat_capa_l%ptr(:,isoil,1) * soil%bulk_dens_l%ptr(:,isoil,1) * soil%soil_depth_l%ptr(:,isoil,1)) 
          ENDIF 
       ENDDO 
  
       !> 2.3 update net radiation flux at the surface, and implied soil temperature
       !!
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

    !> 3.0 update soil moisture
    !!
    SELECT CASE(soil_config%soil_moisture_scheme)
       !> 3.1 bucket soil water scheme
       !!
       CASE("bucket")
          !> 3.1.1 update soil moisture
          !! dS/dt = rain - transpiration - evaporation
          soil%w_soil_root%ptr(1,1) = soil%w_soil_root%ptr(1,1) + (rain + snow &
               - veg%transpiration%ptr - evaporation) / 1000._wp * dtime
   
          !> 3.1.2 calculate drainage if any [mm/s]
          !!
          WHERE(soil%w_soil_root%ptr > soil%w_soil_root_fc%ptr)
            soil%drainage%ptr = ( soil%w_soil_root%ptr - soil%w_soil_root_fc%ptr ) * 1000._wp/ dtime
          ELSEWHERE
            soil%drainage%ptr = 0.0_wp
          ENDWHERE
          soil%w_soil_root%ptr = soil%w_soil_root%ptr - soil%drainage%ptr / 1000.0_wp * dtime
          soil%w_soil_root%ptr = MAX(soil%w_soil_root%ptr,0.0_wp)
   
          !> 3.1.3 calculate plant fractional water availability 
          !!
          soil%w_soil_root_theta%ptr = ( soil%w_soil_root%ptr - soil%w_soil_root_pwp%ptr ) / &
                                       ( soil%w_soil_root_fc%ptr - soil%w_soil_root_pwp%ptr ) 
          soil%w_soil_root_theta%ptr = MAX(soil%w_soil_root_theta%ptr,0.0_wp)

          !> 3.1.4 calculate soil wtaer potential according to Campbell et al 1974 
          soil%w_soil_root_pot%ptr = k_soil_phi1 * & 
              ( soil%w_soil_root%ptr / soil%w_soil_root_fc%ptr )**k_soil_phi2

          !> 3.1.5 update soil moisture by layer 
          !! using mean theta of the rooting zone for those layer that contain roots
          !! note that for ease of computation this is not 100% mass conserving and not meant to be!
          !!
          DO isoil = 1,soil_config%nsoil_water
             WHERE(soil%root_fraction_l%ptr(:,isoil,:)>eps8)
                soil%w_soil_sl%ptr(:,isoil,:) = soil%w_soil_pwp_sl%ptr(:,isoil,:) + & 
                   soil%w_soil_root_theta%ptr * (soil%w_soil_fc_sl%ptr(:,isoil,:) - soil%w_soil_pwp_sl%ptr(:,isoil,:))
             ENDWHERE
          ENDDO
    
       !> 3.2 multi-layer soil moisture scheme 
       !!
       CASE("layer")

          ! save current water levels for percolation calculation later
          soil%percolation_sl%ptr = soil%w_soil_sl%ptr

          !> 3.2.1 determine fraction of evaporation per layer, limited to the first two soil layers
          !!
          fevaporation(:) = 0.0_wp
          IF(soil%w_soil_sl%ptr(1,1,1) > eps8) THEN
             fevaporation(1) = MAX(1.0_wp,evaporation/soil%w_soil_sl%ptr(1,1,1))
          ENDIF
          fevaporation(2) = 1._wp - fevaporation(1) 

          !> 3.2.2 determine fraction of transpiration per layer
          !!
          ftranspiration(:) = 0.0_wp
          DO isoil = 1,soil_config%nsoil_water
             IF(soil%w_soil_fc_sl%ptr(1,isoil,1) - soil%w_soil_pwp_sl%ptr(1,isoil,1) > eps8)THEN
                hlp1 = (soil%w_soil_sl%ptr(1,isoil,1) - soil%w_soil_pwp_sl%ptr(1,isoil,1)) / & 
                       (soil%w_soil_fc_sl%ptr(1,isoil,1) - soil%w_soil_pwp_sl%ptr(1,isoil,1))
                hlp1 = calc_soil_moisture_stress(hlp1,lctlib%lambda_soil_gs,lctlib%weibull_k_soil_gs)
                ftranspiration(isoil) = soil%root_fraction_l%ptr(1,isoil,1) * hlp1 
             ENDIF
          ENDDO
          IF(sum(ftranspiration(:)) > eps8 ) THEN 
             ftranspiration(:) = ftranspiration(:) / sum(ftranspiration(:))
          ELSE 
             ftranspiration(:) = 0.0_wp
          ENDIF

          !> 3.2.3 determine water balance of the top layer 
          !! infiltration rate (limited to attaining field capacity), surface runoff and deep drainage
          !! a fraction of runoff generation is leaked below the first layer
          !! a fraction of infiltration is leaked below the first layer (preferential flow) 
          ! first guess
          infiltration = MAX(MIN(soil%w_soil_fc_sl%ptr(1,1,1) - soil%w_soil_sl%ptr(1,1,1), & 
                         ( rain + snow ) / 1000._wp * dtime),0.0_wp)
          soil%runoff%ptr(:,:) = rain + snow - infiltration * 1000._wp / dtime

          ! modify surface runoff and infiltration such that a constant fraction is actually leaked below the first layer
          soil%drainage_sl%ptr(:,1,:) = drain_frac / soil%soil_depth_l%ptr(:,1,:) * infiltration + fdrain_runoff * soil%runoff%ptr(:,:)
          soil%runoff%ptr(:,:) = ( 1.0_wp - fdrain_runoff ) * soil%runoff%ptr(:,:)  
          infiltration = ( 1.0_wp - drain_frac / soil%soil_depth_l%ptr(1,1,1) ) * infiltration

          soil%w_soil_sl%ptr(:,1,:) = soil%w_soil_sl%ptr(:,1,:) + (infiltration - &
               evaporation * fevaporation(1) + veg%transpiration%ptr * ftranspiration(1)) / 1000._wp * dtime

          !> 3.2.4 update deeper soil layers given infiltration, transpiration and evaporation
          !!
          DO isoil = 2,soil_config%nsoil_water

               ! drainage from this layer due to preferential flow of infiltering water flux
               soil%drainage_sl%ptr(:,isoil,:) = drain_frac / soil%soil_depth_l%ptr(:,isoil,:) * & 
                                                 soil%drainage_sl%ptr(:,isoil-1,:)

               ! water balance of soil layer
               soil%w_soil_sl%ptr(:,isoil,:) = soil%w_soil_sl%ptr(:,isoil,:) + & 
                    (( 1.0_wp - drain_frac / soil%soil_depth_l%ptr(:,isoil,:) ) * soil%drainage_sl%ptr(:,isoil-1,:) - & 
                     evaporation * fevaporation(isoil) + veg%transpiration%ptr * ftranspiration(isoil)) / 1000._wp * dtime

               ! drain surplus water to lower layer to keep water storage below field capacity
               WHERE(soil%w_soil_sl%ptr(:,isoil,:) > soil%w_soil_fc_sl%ptr(:,isoil,:))
                  soil%drainage_sl%ptr(:,isoil,:) = soil%drainage_sl%ptr(:,isoil,:) + & 
                                                    (soil%w_soil_sl%ptr(:,isoil,:) - soil%w_soil_fc_sl%ptr(:,isoil,:))*1000._wp / dtime
                  soil%w_soil_sl%ptr(:,isoil,:) = soil%w_soil_fc_sl%ptr(:,isoil,:) 
               ENDWHERE 

          ENDDO

          !> 3.2.5 update deeper soil layer by diffusion
          !!
          WHERE(soil%w_soil_fc_sl%ptr(1,:,1)>eps8)
             theta_sl(:) = soil%w_soil_sl%ptr(1,:,1) / soil%w_soil_fc_sl%ptr(1,:,1)
          ENDWHERE          
          DO isoil = 2,soil_config%nsoil_water
             ! diffusivity is taken here as linerarily dependent on the fractional water content of a layer
             ! diffusion is taken to be proportional to the fractional water content difference 
             hlp1 = kdiff_theta * theta_sl(isoil) * dtime * & 
                   (theta_sl(isoil-1) - theta_sl(isoil)) / &  
                   ((soil%soil_depth_l%ptr(1,isoil-1,1)+(soil%soil_depth_l%ptr(1,isoil,1)))/2.0_wp)

             soil%w_soil_sl%ptr(:,isoil-1,:) = soil%w_soil_sl%ptr(:,isoil-1,:) - hlp1
             soil%w_soil_sl%ptr(:,isoil,:) = soil%w_soil_sl%ptr(:,isoil,:) + hlp1 
             soil%drainage_sl%ptr(:,isoil-1,:) = soil%drainage_sl%ptr(:,isoil-1,:) + hlp1 * 1000._wp / dtime 

          ENDDO

          !> 3.2.6 calculate rooting zone water content and theta for compatability with bucket scheme and vegetation water stress
          !!
          soil%w_soil_root%ptr = soil%w_soil_sl%ptr(:,1,:) 
          DO isoil=2,soil_config%nsoil_water
             WHERE(soil%root_depth%ptr > soil_config%soil_depth_water(isoil-1))
                WHERE(soil%root_depth%ptr >= soil_config%soil_depth_water(isoil))
                    soil%w_soil_root%ptr = soil%w_soil_root%ptr + soil%w_soil_sl%ptr(:,isoil,:) 
                ELSEWHERE
                    soil%w_soil_root%ptr = soil%w_soil_root%ptr + soil%w_soil_sl%ptr(:,isoil,:) * &  
                       (soil%root_depth%ptr-soil_config%soil_depth_water(isoil-1)) / & 
                       (soil_config%soil_depth_water(isoil)-soil_config%soil_depth_water(isoil-1))
                ENDWHERE
             ENDWHERE
          ENDDO
          soil%w_soil_root_theta%ptr = ( soil%w_soil_root%ptr - soil%w_soil_root_pwp%ptr ) / &
                                       ( soil%w_soil_root_fc%ptr - soil%w_soil_root_pwp%ptr ) 
          soil%w_soil_root_theta%ptr = MAX(soil%w_soil_root_theta%ptr,0.0_wp)

          !> 3.2.7. calculate root zone water potential
          !!
          soil%w_soil_root_pot%ptr = 0.0_wp
          DO isoil=1,soil_config%nsoil_water
             WHERE(soil%w_soil_fc_sl%ptr(:,1,:) > eps8)
                 soil%w_soil_root_pot%ptr = soil%w_soil_root_pot%ptr + soil%root_fraction_l%ptr(:,isoil,:) * &
                          k_soil_phi1 * ( MAX(min_theta_wpot,soil%w_soil_sl%ptr(:,isoil,:) / soil%w_soil_fc_sl%ptr(:,isoil,:)) )**k_soil_phi2
             ENDWHERE
          ENDDO

          !> 3.2.8. diagnose deep drainage
          !!
          soil%drainage%ptr = soil%drainage_sl%ptr(:,soil_config%nsoil_water,:) 
          WHERE(soil%percolation_sl%ptr > eps8)
             soil%percolation_sl%ptr = MIN(MAX((soil%drainage_sl%ptr/1000_wp * dtime )/ soil%percolation_sl%ptr,-1.0_wp),1.0_wp)*1000._wp/dtime
          ELSEWHERE
             soil%percolation_sl%ptr = 0.0_wp
          ENDWHERE

    ENDSELECT

    veg%dphi%ptr = soil%w_soil_root_pot%ptr(1,1) - lctlib%phi_leaf_min - grav * rhoh2o * veg%height%ptr * 1.e-6_wp

  END SUBROUTINE update_soil_physics



END MODULE mo_soil_update_physics
