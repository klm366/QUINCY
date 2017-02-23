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
  SUBROUTINE update_soil_physics(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_base_constants,        ONLY: dtime
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,      ONLY: r_gas_dryair, LatentHeatVaporization, stbo, zemiss_def
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
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)    :: options          !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)       :: &                    ! 1D (ics-ice)
                                             q_sat        , &
                                             q_sat_soil   , &
                                             evaporation  , &
                                             infiltration , &
                                             hlp1
    REAL(wp),DIMENSION(options%nc,soil_config%nsoil_water) :: &   ! 2D (ics-ice,z-axis)
                                             ftranspiration , &
                                             fevaporation   , &
                                             theta_sl
    INTEGER                               :: isoil
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_physics'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)                   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(LAND2ATM_)                   !     TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(SOIL_)                       !     TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk                :: t_air                  ! ATM2LAND_
    dsl4jsb_Real2D_onChunk                :: q_air
    dsl4jsb_Real2D_onChunk                :: press_srf
    dsl4jsb_Real2D_onChunk                :: ga
    dsl4jsb_Real2D_onChunk                :: lw_srf_down
    dsl4jsb_Real2D_onChunk                :: rain
    dsl4jsb_Real2D_onChunk                :: snow
    dsl4jsb_Real2D_onChunk                :: drag_srf
    dsl4jsb_Real2D_onChunk                :: net_radiation          ! LAND2ATM_
    dsl4jsb_Real2D_onChunk                :: ground_heat_flx
    dsl4jsb_Real2D_onChunk                :: latent_heat_flx
    dsl4jsb_Real2D_onChunk                :: sensible_heat_flx
    dsl4jsb_Real2D_onChunk                :: root_depth             ! SOIL_
    dsl4jsb_Real2D_onChunk                :: w_soil_root
    dsl4jsb_Real2D_onChunk                :: w_soil_root_pwp
    dsl4jsb_Real2D_onChunk                :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk                :: w_soil_root_theta
    dsl4jsb_Real2D_onChunk                :: w_soil_root_pot
    dsl4jsb_Real2D_onChunk                :: runoff
    dsl4jsb_Real2D_onChunk                :: drainage
    dsl4jsb_Real3D_onChunk                :: soil_depth_l
    dsl4jsb_Real3D_onChunk                :: bulk_dens_l
    dsl4jsb_Real3D_onChunk                :: heat_capa_l
    dsl4jsb_Real3D_onChunk                :: therm_cond_l
    dsl4jsb_Real3D_onChunk                :: root_fraction_l
    dsl4jsb_Real3D_onChunk                :: w_soil_pwp_sl
    dsl4jsb_Real3D_onChunk                :: w_soil_fc_sl 
    dsl4jsb_Real3D_onChunk                :: w_soil_sl
    dsl4jsb_Real3D_onChunk                :: percolation_sl
    dsl4jsb_Real3D_onChunk                :: drainage_sl
    dsl4jsb_Real3D_onChunk                :: t_soil
    dsl4jsb_Real2D_onChunk                :: canopy_cond            ! ASSIMI_
    dsl4jsb_Real2D_onChunk                :: sw_srf_net             ! RAD_
    dsl4jsb_Real2D_onChunk                :: transpiration          ! VEG_
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(LAND2ATM_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_air)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_air)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, press_srf)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, ga)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, lw_srf_down)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, rain)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, snow)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, drag_srf)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(LAND2ATM_, net_radiation)
    dsl4jsb_Get_var2D_onChunk(LAND2ATM_, ground_heat_flx)
    dsl4jsb_Get_var2D_onChunk(LAND2ATM_, latent_heat_flx)
    dsl4jsb_Get_var2D_onChunk(LAND2ATM_, sensible_heat_flx)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, root_depth)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_pwp)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_fc)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_theta)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_pot)
    dsl4jsb_Get_var2D_onChunk(SOIL_, runoff)
    dsl4jsb_Get_var2D_onChunk(SOIL_, drainage)
    dsl4jsb_Get_var3D_onChunk(SOIL_, soil_depth_l)
    dsl4jsb_Get_var3D_onChunk(SOIL_, bulk_dens_l)
    dsl4jsb_Get_var3D_onChunk(SOIL_, heat_capa_l)
    dsl4jsb_Get_var3D_onChunk(SOIL_, therm_cond_l)
    dsl4jsb_Get_var3D_onChunk(SOIL_, root_fraction_l)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_pwp_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_fc_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, percolation_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, drainage_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, t_soil)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(RAD_, sw_srf_net)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, transpiration)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Latent heat flux
    !!

    !> 1.1 calculate saturating specific humidity 
    !!
    q_sat(:) = calc_spec_humidity_sat(t_air(:),press_srf(:))
    q_sat_soil(:) = calc_spec_humidity_sat(t_soil(:,1),press_srf(:))

    !> 1.2. calulate transpiration (kg / m2 / s)
    !! T = rho * (qsat - qair ) * Gtot in mm/s 
    !!
    WHERE (canopy_cond(:) > eps8)
       transpiration(:) = press_srf(:) / ( r_gas_dryair * t_air(:) ) * & 
                  (q_sat(:) - q_air(:)) / (1._wp / canopy_cond(1) + 1._wp / ga(:) )  
    ELSE WHERE
       transpiration(:) = 0.0_wp 
    END WHERE

    !> 1.3 soil evaporation (kg / m2 / s)
    !! E = rho * (qsat - qair) * theta * 1/ga_soil (taken as ga_soil(:) = ga(:) / 10)
    !!
    SELECT CASE(soil_config%soil_moisture_scheme)
       CASE("bucket") ! ! Case of bucket water scheme
          evaporation(:) = press_srf(:) / ( r_gas_dryair * t_air(:) ) * & 
                        (q_sat(:) - q_air(:)) / ( 100._wp / ga(:)) * & 
                        w_soil_root_theta(:)
          WHERE(sw_srf_net(:) < eps8) evaporation(:) = 0.0_wp
       
       CASE("layer") ! Case of layer scheme
          evaporation(:) = press_srf(:) / ( r_gas_dryair * t_air(:) ) * & 
                        (q_sat - q_air(:)) / ( 10._wp / ga(:)) * & 
                        SUM(w_soil_sl(:,1:2)) / SUM(w_soil_fc_sl(:,1:2)) 
          WHERE(evaporation(:) < eps8) evaporation(:) = 0.0_wp
    ENDSELECT

    !> 1.4 calculate latent heat flux (W m-2) 
    !!
    latent_heat_flx(:) = LatentHeatVaporization * (evaporation(:) + transpiration(:))

    !> 2.0 Estimate surface temperature and sensible heat flux 
    !!
    IF (soil_config%flag_calc_soil_energy) THEN ! case of a simple diagnostic soil energy calculation
       !> 2.1 simple estimate of sensible heat flux
       !! (ignoring instabilities and actual surface moisture...)
       !!
       sensible_heat_flx(:) = heat_transfer_coef(drag_srf(:), model%options%steplen) * & 
                             (surface_dry_static_energy(t_soil(:,1),q_sat_soil(:)) - & 
                              surface_dry_static_energy(t_air(:),q_air(:)))

       !> 2.2 Estimate ground flux and update deep soil temperature
       !! working upwards to estimate net transfer implied by previous temperature profile 
       !! (ignores interaction of diffusion with water or ice)
       !!
       DO isoil = soil_config%nsoil_water,2,-1
          WHERE(soil_depth_l(:,isoil) > eps8)
             ground_heat_flx(:) = (t_soil(:,isoil-1)-t_soil(:,isoil)) * therm_cond_l(:,isoil) / & 
                        (soil_depth_l(:,isoil)/2.0_wp+soil_depth_l(:,isoil-1)/2.0_wp)
             t_soil(:,isoil) = t_soil(:,isoil) + ground_heat_flx(:) * dtime / & 
                        (heat_capa_l(:,isoil) * bulk_dens_l(:,isoil) * soil_depth_l(:,isoil)) 
          END WHERE 
       ENDDO 
  
       !> 2.3 update net radiation flux at the surface, and implied soil temperature
       !!
       net_radiation(:) = latent_heat_flx(:) + sensible_heat_flx(:) - ground_heat_flx(:) 
       t_soil(:,1) = t_soil(:,1) + ( sw_srf_net(:) - net_radiation(:) + & 
                                ! simple estimate of surface net longwave flux
                                lw_srf_down(:) - stbo * zemiss_def * t_soil(:,1)**4 ) / & 
                                (heat_capa_l(:,1) * bulk_dens_l(:,1) * soil_depth_l(:,1)) * &
                                dtime 

    ELSE ! diagnose sensible heat flux from net radiation and set soil temperature to air temperature 
       ground_heat_flx(:) = 0.0_wp 
       sensible_heat_flx(:) = net_radiation(:) - latent_heat_flx(:) 
       t_soil(:,1) = t_air(:)
    ENDIF

    !> 3.0 update soil moisture
    !!
    SELECT CASE(soil_config%soil_moisture_scheme)
       !> 3.1 bucket soil water scheme
       !!
       CASE("bucket")
          !> 3.1.1 update soil moisture
          !! dS/dt = rain - transpiration - evaporation
          w_soil_root(:) = w_soil_root(:) + (rain(:) + snow(:) &
               - transpiration(:) - evaporation(:)) / 1000._wp * dtime
   
          !> 3.1.2 calculate drainage if any [mm/s]
          !!
          WHERE(w_soil_root(:) > w_soil_root_fc(:))
            drainage(:) = ( w_soil_root(:) - w_soil_root_fc(:) ) * 1000._wp/ dtime
          ELSEWHERE
            drainage(:) = 0.0_wp
          END WHERE
          w_soil_root(:) = w_soil_root(:) - drainage(:) / 1000.0_wp * dtime
          w_soil_root(:) = MAX(w_soil_root(:),0.0_wp)
   
          !> 3.1.3 calculate plant fractional water availability 
          !!
          w_soil_root_theta(:) = ( w_soil_root(:) - w_soil_root_pwp(:) ) / &
                                       ( w_soil_root_fc(:) - w_soil_root_pwp(:) ) 
          w_soil_root_theta(:) = MAX(w_soil_root_theta(:),0.0_wp)

          !> 3.1.4 calculate soil wtaer potential according to Campbell et al 1974 
          w_soil_root_pot(:) = k_soil_phi1 * & 
              ( w_soil_root(:) / w_soil_root_fc(:) )**k_soil_phi2

          !> 3.1.5 update soil moisture by layer 
          !! using mean theta of the rooting zone for those layer that contain roots
          !! note that for ease of computation this is not 100% mass conserving and not meant to be!
          !!
          DO isoil = 1,soil_config%nsoil_water
             WHERE(root_fraction_l(:,isoil) > eps8)
                w_soil_sl(:,isoil) = w_soil_pwp_sl(:,isoil) + & 
                   w_soil_root_theta(:) * (w_soil_fc_sl(:,isoil) - w_soil_pwp_sl(:,isoil))
             END WHERE
          ENDDO
    
       !> 3.2 multi-layer soil moisture scheme 
       !!
       CASE("layer")

          ! save current water levels for percolation calculation later
          percolation_sl(:,:) = w_soil_sl(:,:)

          !> 3.2.1 determine fraction of evaporation per layer, limited to the first two soil layers
          !!
          fevaporation(:,:) = 0.0_wp
          WHERE(w_soil_sl(:,1) > eps8)
             fevaporation(:,1) = MAX(1.0_wp,evaporation(:)/w_soil_sl(:,1))
          END WHERE

          !> 3.2.2 determine fraction of transpiration per layer
          !!
          ftranspiration(:,:) = 0.0_wp
          DO isoil = 1,soil_config%nsoil_water
             WHERE(w_soil_fc_sl(:,isoil) - w_soil_pwp_sl(:,isoil) > eps8)
                hlp1(:) = (w_soil_sl(:,isoil) - w_soil_pwp_sl(:,isoil)) / & 
                       (w_soil_fc_sl(:,isoil) - w_soil_pwp_sl(:,isoil))
                hlp1(:) = calc_soil_moisture_stress(hlp1(:),lctlib%lambda_soil_gs,lctlib%weibull_k_soil_gs)
                ftranspiration(:,isoil) = root_fraction_l(:,isoil) * hlp1(:) 
             END WHERE
          ENDDO
          IF(sum(ftranspiration(:,:)) > eps8 ) THEN 
             ftranspiration(:,:) = ftranspiration(:,:) / sum(ftranspiration(:,:))
          ELSE 
             ftranspiration(:,:) = 0.0_wp
          ENDIF

          !> 3.2.3 determine water balance of the top layer 
          !! infiltration rate (limited to attaining field capacity), surface runoff and deep drainage
          !! a fraction of runoff generation is leaked below the first layer
          !! a fraction of infiltration is leaked below the first layer (preferential flow) 
          ! first guess
          infiltration(:) = MAX(MIN(w_soil_fc_sl(:,1) - w_soil_sl(:,1), & 
                         ( rain(:) + snow(:) ) / 1000._wp * dtime),0.0_wp)
          runoff(:) = rain(:) + snow(:) - infiltration(:) * 1000._wp / dtime

          ! modify surface runoff and infiltration such that a constant fraction is actually leaked below the first layer
          drainage_sl(:,1) = drain_frac / soil_depth_l(:,1) * infiltration(:) + fdrain_runoff * runoff(:)
          runoff(:) = ( 1.0_wp - fdrain_runoff ) * runoff(:)  
          infiltration(:) = ( 1.0_wp - drain_frac / soil_depth_l(:,1) ) * infiltration(:)

          w_soil_sl(:,1) = w_soil_sl(:,1) + (infiltration(:) - &
               evaporation(:) * fevaporation(:,1) + transpiration(:) * ftranspiration(:,1)) / 1000._wp * dtime

          WHERE(w_soil_sl(:,1) > w_soil_fc_sl(:,1))
               runoff(:) = runoff(:) + (w_soil_sl(:,1) - w_soil_fc_sl(:,1))*1000._wp / dtime
               w_soil_sl(:,1) = w_soil_fc_sl(:,1) 
          END WHERE 


          !> 3.2.4 update deeper soil layers given infiltration, transpiration and evaporation
          !!
          DO isoil = 2,soil_config%nsoil_water

               ! drainage from this layer due to preferential flow of infiltering water flux
               drainage_sl(:,isoil) = drain_frac / soil_depth_l(:,isoil) * & 
                                                 drainage_sl(:,isoil-1)

               ! water balance of soil layer
               w_soil_sl(:,isoil) = w_soil_sl(:,isoil) + & 
                    (( 1.0_wp - drain_frac / soil_depth_l(:,isoil) ) * drainage_sl(:,isoil-1) - & 
                     evaporation(:) * fevaporation(:,isoil) + transpiration(:) * ftranspiration(:,isoil)) / 1000._wp * dtime

               ! drain surplus water to lower layer to keep water storage below field capacity
               WHERE(w_soil_sl(:,isoil) > w_soil_fc_sl(:,isoil))
                  drainage_sl(:,isoil) = drainage_sl(:,isoil) + & 
                                                    (w_soil_sl(:,isoil) - w_soil_fc_sl(:,isoil))*1000._wp / dtime
                  w_soil_sl(:,isoil) = w_soil_fc_sl(:,isoil) 
               END WHERE 

          ENDDO

          !> 3.2.5 update deeper soil layer by diffusion
          !!
          WHERE(w_soil_fc_sl(:,:) > eps8)
             theta_sl(:,:) = w_soil_sl(:,:) / w_soil_fc_sl(:,:)
          END WHERE          
          DO isoil = 2,soil_config%nsoil_water
             ! diffusivity is taken here as linerarily dependent on the fractional water content of a layer
             ! diffusion is taken to be proportional to the fractional water content difference 
             hlp1 = kdiff_theta * theta_sl(:,isoil) * dtime * & 
                   (theta_sl(:,isoil-1) - theta_sl(:,isoil)) / &  
                   ((soil_depth_l(:,isoil-1)+(soil_depth_l(:,isoil)))/2.0_wp)

             w_soil_sl(:,isoil-1) = w_soil_sl(:,isoil-1) - hlp1(:)
             w_soil_sl(:,isoil) = w_soil_sl(:,isoil) + hlp1(:) 
             drainage_sl(:,isoil-1) = drainage_sl(:,isoil-1) + hlp1(:) * 1000._wp / dtime 

          ENDDO

          !> 3.2.6 calculate rooting zone water content and theta for compatability with bucket scheme and vegetation water stress
          !!
          w_soil_root(:) = w_soil_sl(:,1) 
          DO isoil=2,soil_config%nsoil_water
             WHERE(root_depth(:) > soil_config%soil_depth_water(isoil-1))
                WHERE(root_depth(:) >= soil_config%soil_depth_water(isoil))
                    w_soil_root(:) = w_soil_root(:) + w_soil_sl(:,isoil) 
                ELSEWHERE
                    w_soil_root(:) = w_soil_root(:) + w_soil_sl(:,isoil) * &  
                       (root_depth(:)-soil_config%soil_depth_water(isoil-1)) / & 
                       (soil_config%soil_depth_water(isoil)-soil_config%soil_depth_water(isoil-1))
                ENDWHERE
             ENDWHERE
          ENDDO
          w_soil_root_theta(:) = ( w_soil_root(:) - w_soil_root_pwp(:) ) / &
                                       ( w_soil_root_fc(:) - w_soil_root_pwp(:) ) 
          w_soil_root_theta(:) = MAX(w_soil_root_theta(:),0.0_wp)

          !> 3.2.7. calculate root zone water potential
          !!
          w_soil_root_pot(:) = 0.0_wp
          DO isoil=1,soil_config%nsoil_water
             WHERE(w_soil_fc_sl(:,1) > eps8)
                 w_soil_root_pot(:) = w_soil_root_pot(:) + root_fraction_l(:,isoil) * &
                          k_soil_phi1 * ( MAX(min_theta_wpot,w_soil_sl(:,isoil) / w_soil_fc_sl(:,isoil)) )**k_soil_phi2
             ENDWHERE
          ENDDO

          !> 3.2.8. diagnose deep drainage
          !!
          drainage(:) = drainage_sl(:,soil_config%nsoil_water) 
          WHERE(percolation_sl(:,:) > eps8)
             percolation_sl(:,:) = MIN(MAX((drainage_sl(:,:)/1000_wp * dtime )/ percolation_sl(:,:),-1.0_wp),1.0_wp)*1000._wp/dtime
          ELSEWHERE
             percolation_sl(:,:) = 0.0_wp
          ENDWHERE

    ENDSELECT


  END SUBROUTINE update_soil_physics



END MODULE mo_soil_update_physics
