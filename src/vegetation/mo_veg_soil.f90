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
         plant%f_n_demand = MAX(plant%moving_average%growth_req_n_week - nlab / clab, 0.0_wp) 
      ELSE
         plant%f_n_demand = 1.0_wp
      ENDIF
      IF(nlab > eps8) THEN
         plant%f_p_demand = MAX(plant%moving_average%growth_req_p_week - plab / nlab, 0.0_wp)
      ELSE
         plant%f_p_demand = 1.0_wp
      ENDIF

      !> 1.3 calculate N, P uptake according to OCN 
      !! 
      !!
      plant%flux%uptake_n = lctlib%vmax_uptake_n * &
              plant%moving_average%demand_uptake_n * & 
              soil%ammonium * ( lctlib%soil_n_aff1 + & 
                1._wp / ( lctlib%soil_n_aff2 + soil%ammonium )) * & 
              plant%pool%fine_root%carbon
      plant%flux%uptake_p = lctlib%vmax_uptake_p * & 
              plant%moving_average%demand_uptake_p * & 
              soil%phosphate * ( lctlib%soil_p_aff1 + & 
                1._wp / ( lctlib%soil_p_aff2 + soil%phosphate )) * & 
              plant%pool%fine_root%carbon

    !> 2.0 ELSE: maintain prescribed C:N:P irrespective of nutrient uptake capacity == no nutrient limitation case 
    !! 
    ELSE
      ! uptake required to maintain a presribed C:N ratio of the labile pool in µmol N / m2 / s
      plant%flux%uptake_n = MAX(plant%pool%labile%carbon / lctlib%cn_leaf - & 
                                        plant%pool%labile%nitrogen, 0.0_wp) * 1.e6 / dt
      ! uptake required to maintain a presribed N:P ratio of the labile pool in µmol P / m2 / s
      plant%flux%uptake_p = MAX(plant%pool%labile%nitrogen / lctlib%np_leaf + & 
                                           plant%flux%uptake_n * dt / 1.e6 - &
                                           plant%pool%labile%phosphorus, 0.0_wp) * 1.e6 / dt

    ENDIF     

    !> 3.0 Calculate potential N uptake per unit root
    !!
    plant%flux%unit_uptake_n = lctlib%vmax_uptake_n * &
                               plant%moving_average%demand_uptake_n * & 
                               soil%ammonium * ( lctlib%soil_n_aff1 + & 
                               1._wp / ( lctlib%soil_n_aff2 + soil%ammonium ))
    plant%flux%unit_uptake_p = lctlib%vmax_uptake_p * &
                               plant%moving_average%demand_uptake_p * & 
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
    soil%flux%loss%litter = 1.0_wp/tau_litter * soil%pool%litter
    
    !> 2 Update soil litter pool. Probably needs splitting into structural and non-structural litter
    !!
    soil%pool%litter = soil%pool%litter + plant%flux%litterfall%leaf
    soil%pool%litter = soil%pool%litter + plant%flux%litterfall%fine_root
    soil%pool%litter = soil%pool%litter + plant%flux%litterfall%coarse_root
    soil%pool%litter = soil%pool%litter + plant%flux%litterfall%heart_wood

    soil%pool%litter = soil%pool%litter - soil%flux%loss%litter * dt
    
    !> 3 Update mineral pool
    !! 
!SZ BLARPP. check units
    soil%ammonium = soil%ammonium   + soil%flux%loss%litter%nitrogen*dt*  1.e3 / soil_depth * 1.e-3
    soil%ammonium = soil%ammonium   - plant%flux%uptake_n * dt/1.e6 * 1.e3 / soil_depth * 1.e-3 
    soil%ammonium = soil%ammonium   + lctlib%bnf_base * dt * 1.e3 / soil_depth * 1.e-3
    !soil%phosphate = soil%phosphate + soil%flux%loss%litter%phosphorus*dt*  1.e3 / soil_depth * 1.e-3
    !soil%phosphate = soil%phosphate - plant%flux%uptake_p * dt/1.e6 * 1.e3 / soil_depth * 1.e-3 

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

    USE mo_base_memory,       ONLY: t_slm_soil, t_slm_tile
    USE mo_veg_memory,        ONLY: t_slm_plant
    USE mo_base_constants,    ONLY: dt
    USE mo_atmland_memory,    ONLY: atm2land
    USE mo_atmland_util,      ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,  ONLY: r_gas_dryair

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile     !< one tile
    ! 0.2 Local
    TYPE(t_slm_soil),POINTER                :: soil     ! soil variables per tile
    TYPE(t_slm_plant),POINTER               :: plant    ! plant variables per tile
    INTEGER  :: t
    REAL(wp) :: q_sat
    REAL(wp) :: evaporation

    !  bi-directional pointer
    plant => tile%plant
    soil  => tile%soil

    IF(veg_config%flag_calc_soil_moisture) THEN ! case of dynamic water uptake calculation 

       ! saturating specific humidity 

       q_sat = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf)

       IF(plant%gc > eps8)THEN
          ! T = rho * (qsat - qair ) * Gtot in mm/s 
          plant%flux%transpiration = atm2land%press_srf / ( r_gas_dryair * atm2land%t_air ) * & 
                     (q_sat - atm2land%q_air) / (1._wp / plant%gc + 1._wp / atm2land%ga )  
       ELSE 
          plant%flux%transpiration = 0.0_wp 
       ENDIF
       ! E = rho * (qsat - qair) * theta * 1/ga_soil (taken as ga_soil = ga / 10)
       evaporation = atm2land%press_srf / ( r_gas_dryair * atm2land%t_air ) * & 
                     (q_sat - atm2land%q_air) * soil%theta / 100.0_wp / ( 10._wp / atm2land%ga)

       ! update soil moisture dS = rain - transpiration - evaporation 
       soil%moisture = soil%moisture + (atm2land%rain + atm2land%snow & 
            - plant%flux%transpiration - evaporation) * dt
       ! constrain to possible solutions
       soil%moisture = MAX(MIN(soil%moisture,veg_config%soil_awc_prescribe),0.0_wp)
       ! scale to relative water availability (%)
       soil%theta = soil%moisture / veg_config%soil_awc_prescribe * 100._wp

   ELSE 
       ! use prescribed value of theta from namelist
       soil%theta = veg_config%soil_theta_prescribe 
 
       ! calculate corresponding soil moisture value (mm)     
       soil%moisture = veg_config%soil_awc_prescribe * soil%theta / 100.0_wp

   ENDIF 
  END SUBROUTINE update_vegsoil_physics

END MODULE mo_veg_soil
