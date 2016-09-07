!>
!! @brief diverse subroutines and functions to mimick soil for the stand-alone plant-model
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-03-31)
!!
MODULE mo_veg_soil

  USE mo_kind,               ONLY: wp, eps8
  USE mo_jsb_lctlib_class,   ONLY: t_lctlib_element
  USE mo_veg_config_class,   ONLY: veg_config
  USE mo_veg_constants


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_vegsoil_bgc, update_vegsoil_physics, update_soil_pool


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> routine to simulate N and P uptake from soil in the absence of a proper soil model
  !! 
  !! @par alternative use
  !! supply N and P to maintain a certain target C:N:P of the labile pool
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_bgc(tile, lctlib)

    USE mo_base_memory,      ONLY: t_slm_soil, t_slm_tile
    USE mo_veg_memory,       ONLY: t_slm_plant
    USE mo_base_constants,   ONLY: dt

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)        :: lctlib        !< land-cover-type library - parameters for one pft
    TYPE(t_slm_tile),INTENT(inout), TARGET   :: tile          !< one tile
    ! 0.2 Local
    TYPE(t_slm_soil),POINTER      :: soil               ! soil variables per tile
    TYPE(t_slm_plant),POINTER     :: plant              ! plant variables per tile
    REAL(wp) :: clab, nlab, plab                        ! plant labile pools of carbon, nitrogen and phosphorus

    !  bi-directional pointer
    plant => tile%plant
    soil  => tile%soil


    !> 1.0 IF: root and demand based N/P uptake   
    !! 
    !!
    IF(veg_config%flag_calc_soil_uptake) THEN ! case of root and demand based N/P uptake

      !> 1.1 calculate current "labile" C:N:P   
      !! 
      clab = plant%pool%labile%carbon
      nlab = plant%pool%labile%nitrogen
      plab = plant%pool%labile%phosphorus

      !> 1.2 convert this into an instantaneous scalar of nutrient demand
      !!     taken as the weekly growth requirement of the nutrient
      !! 
      IF(clab > eps8) THEN
         plant%f_n_demand = MAX(plant%growth_req_n_mavg - nlab / clab, 0.0_wp) 
      ELSE
         plant%f_n_demand = 1.0_wp
      ENDIF
      IF(nlab > eps8) THEN
         plant%f_p_demand = MAX(plant%growth_req_p_mavg - plab / nlab, 0.0_wp)
      ELSE
         plant%f_p_demand = 1.0_wp
      ENDIF

      !> 1.3 calculate N, P uptake according to OCN 
      !! 
      !!
      plant%uptake_n = lctlib%vmax_uptake_n * &
              plant%demand_uptake_n_mavg * & 
              soil%ammonium * ( lctlib%soil_n_aff1 + & 
                1._wp / ( lctlib%soil_n_aff2 + soil%ammonium )) * & 
              plant%pool%fine_root%carbon
      plant%uptake_p = lctlib%vmax_uptake_p * & 
              plant%demand_uptake_p_mavg * & 
              soil%phosphate * ( lctlib%soil_p_aff1 + & 
                1._wp / ( lctlib%soil_p_aff2 + soil%phosphate )) * & 
              plant%pool%fine_root%carbon

    !> 2.0 ELSE: maintain prescribed C:N:P irrespective of nutrient uptake capacity == no nutrient limitation case 
    !! 
    ELSE
      ! uptake required to maintain a presribed C:N ratio of the labile pool in µmol N / m2 / s
      plant%uptake_n = MAX(plant%pool%labile%carbon / lctlib%cn_leaf - & 
                                        plant%pool%labile%nitrogen, 0.0_wp) * 1.e6 / dt
      ! uptake required to maintain a presribed N:P ratio of the labile pool in µmol P / m2 / s
      plant%uptake_p = MAX(plant%pool%labile%nitrogen / lctlib%np_leaf + & 
                                           plant%uptake_n * dt / 1.e6 - &
                                           plant%pool%labile%phosphorus, 0.0_wp) * 1.e6 / dt

    ENDIF     

    !> 3.0 Calculate potential N uptake per unit root
    !!
    plant%unit_uptake_n = lctlib%vmax_uptake_n * &
                               plant%demand_uptake_n_mavg * & 
                               soil%ammonium * ( lctlib%soil_n_aff1 + & 
                               1._wp / ( lctlib%soil_n_aff2 + soil%ammonium ))
    plant%unit_uptake_p = lctlib%vmax_uptake_p * &
                               plant%demand_uptake_p_mavg * & 
                               soil%phosphate * ( lctlib%soil_p_aff1 + & 
                               1._wp / ( lctlib%soil_p_aff2 + soil%phosphate ))

  END SUBROUTINE update_vegsoil_bgc

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> routine to update soil litter and mineral pools
  !! 
  !! @todo add C and P pools, add proper decomposition: mo_veg_soil:update_soil_pool
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_soil_pool(tile, lctlib)

    USE mo_base_memory,      ONLY: t_slm_soil, t_slm_tile
    USE mo_base_constants,   ONLY: dt
    USE mo_veg_memory,       ONLY: t_slm_plant

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)         :: lctlib        !< land-cover-type library - parameters for one pft
    TYPE(t_slm_tile),INTENT(inout), TARGET    :: tile          !< one tile
    ! 0.2 Local
    TYPE(t_slm_soil),POINTER      :: soil               ! soil variables per tile
    TYPE(t_slm_plant),POINTER     :: plant              ! plant variables per tile
    REAl(wp) :: soil_depth = 1.0 ! soil depth, m

    !  bi-directional pointer
    plant => tile%plant
    soil  => tile%soil


    !> 1 Calculate soil litter pool decomposition. Currently only for N
    !!    
    !! Needs temperature and moisture correction
    soil%loss%litter = 1.0_wp/tau_litter * soil%pool%litter
    
    !> 2 Update soil litter pool. Probably needs splitting into structural and non-structural litter
    !!
    soil%pool%litter = soil%pool%litter + plant%litterfall%leaf
    soil%pool%litter = soil%pool%litter + plant%litterfall%fine_root
    soil%pool%litter = soil%pool%litter + plant%litterfall%coarse_root
    soil%pool%litter = soil%pool%litter + plant%litterfall%heart_wood

    soil%pool%litter = soil%pool%litter - soil%loss%litter * dt
    
    !> 3 Update mineral pool
    !! 
!SZ BLARPP. check units
    soil%ammonium = soil%ammonium   + soil%loss%litter%nitrogen*dt*  1.e3 / soil_depth * 1.e-3
    soil%ammonium = soil%ammonium   - plant%uptake_n * dt/1.e6 * 1.e3 / soil_depth * 1.e-3 
    soil%ammonium = soil%ammonium   + lctlib%bnf_base * dt * 1.e3 / soil_depth * 1.e-3
    !soil%phosphate = soil%phosphate + soil%loss%litter%phosphorus*dt*  1.e3 / soil_depth * 1.e-3
    !soil%phosphate = soil%phosphate - plant%uptake_p * dt/1.e6 * 1.e3 / soil_depth * 1.e-3 

    END SUBROUTINE update_soil_pool

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil moisture and theta for the plant only model
  !! 
  !! two cases considered: \n a) dynamic water uptake calculation, \n b) prescribed values from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_physics(tile)

    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
!    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_ --> soil creates conflict

    USE mo_base_memory,           ONLY: t_slm_soil, t_slm_tile
    USE mo_veg_memory,            ONLY: t_slm_plant
    USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    USE mo_land2atm_memory_class, ONLY: t_land2atm_memory

    USE mo_base_constants,        ONLY: dt
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,      ONLY: r_gas_dryair,LatentHeatVaporization

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile     !< one tile
    ! 0.2 Local
    TYPE(t_slm_soil),POINTER                :: soil     ! soil variables per tile
    TYPE(t_slm_plant),POINTER               :: plant    ! plant variables per tile
    INTEGER  :: t
    REAL(wp) :: q_sat
    REAL(wp) :: evaporation
    ! 0.3 Memory allocation
    TYPE(t_atm2land_memory), POINTER :: atm2land__mem   ! dsl4jsb_Def_memory(ATM2LAND_)
    TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    REAL(wp), POINTER                :: t_air
    REAL(wp), POINTER                :: press_srf
    REAL(wp), POINTER                :: q_air
    REAL(wp), POINTER                :: ga
    REAL(wp), POINTER                :: rain
    REAL(wp), POINTER                :: snow
    REAL(wp), POINTER                :: net_radiation
    REAL(wp), POINTER                :: latent_heat_flx
    REAL(wp), POINTER                :: sensible_heat_flx
    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p)            ! dsl4jsb_Get_memory(ATM2LAND_)
      TYPE IS(t_atm2land_memory)
        atm2land__mem => mm
    END SELECT
    SELECT TYPE (mm => tile%mem(LAND2ATM_)%p)
      TYPE IS(t_land2atm_memory)
        land2atm__mem => mm
    END SELECT
    t_air       => atm2land__mem%t_air%ptr
    press_srf   => atm2land__mem%press_srf%ptr
    q_air       => atm2land__mem%q_air%ptr
    ga          => atm2land__mem%ga%ptr
    rain        => atm2land__mem%rain%ptr
    snow        => atm2land__mem%snow%ptr
    net_radiation      => land2atm__mem%net_radiation%ptr
    latent_heat_flx    => land2atm__mem%latent_heat_flx%ptr
    sensible_heat_flx  => land2atm__mem%sensible_heat_flx%ptr


    !  bi-directional pointer
    plant => tile%plant
    soil  => tile%soil

    IF(veg_config%flag_calc_soil_moisture) THEN ! case of dynamic water uptake calculation 

       ! saturating specific humidity 

       q_sat = calc_spec_humidity_sat(t_air,press_srf)

       IF(plant%gc > eps8)THEN
          ! T = rho * (qsat - qair ) * Gtot in mm/s 
          plant%transpiration = press_srf / ( r_gas_dryair * t_air ) * & 
                     (q_sat - q_air) / (1._wp / plant%gc + 1._wp / ga )  
       ELSE 
          plant%transpiration = 0.0_wp 
       ENDIF
       ! E = rho * (qsat - qair) * theta * 1/ga_soil (taken as ga_soil = ga / 10)
       evaporation = press_srf / ( r_gas_dryair * t_air ) * & 
                     (q_sat - q_air) * soil%theta / 100.0_wp / ( 100._wp / ga)
       IF(net_radiation < eps8) evaporation = 0.0_wp

       ! update soil moisture dS = rain - transpiration - evaporation 
       soil%moisture = soil%moisture + (rain + snow & 
            - plant%transpiration - evaporation) * dt
       ! constrain to possible solutions
       soil%moisture = MAX(MIN(soil%moisture,veg_config%soil_awc_prescribe),0.0_wp)
       ! scale to relative water availability (%)
       soil%theta = soil%moisture / veg_config%soil_awc_prescribe * 100._wp

       ! calculate latent heat flux (W m-2) 
       latent_heat_flx = LatentHeatVaporization * (evaporation + plant%transpiration)
       ! diagnose sensible heat flux (shortwave compontent only)
       sensible_heat_flx = net_radiation - latent_heat_flx 

   ELSE 
       ! use prescribed value of theta from namelist
       soil%theta = veg_config%soil_theta_prescribe 
 
       ! calculate corresponding soil moisture value (mm)     
       soil%moisture = veg_config%soil_awc_prescribe * soil%theta / 100.0_wp

   ENDIF 
  END SUBROUTINE update_vegsoil_physics

END MODULE mo_veg_soil
