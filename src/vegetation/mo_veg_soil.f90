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
  USE mo_veg_config_class,   ONLY: veg_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_vegsoil_bgc, update_vegsoil_physics

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> routine to simulate N and P uptake from soil in the absence of a proper soil model
  !! 
  !! @par alternative use
  !! supply N and P to maintain a certain target C:N:P of the labile pool
  !! 
  !! @todo a proper soil model
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_bgc(plant, soil, veg_param)

    USE mo_base_memory,      ONLY: t_slm_soil
    USE mo_veg_memory,       ONLY: t_slm_plant
    USE mo_base_constants,   ONLY: dt
    USE mo_veg_parameters

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_param), INTENT(in) :: veg_param      !< vegetation type parameters
    TYPE(t_slm_soil), INTENT(in)      :: soil           !< soil variables per tile
    TYPE(t_slm_plant), INTENT(inout)  :: plant          !< plant variables per tile
    ! 0.2 Local
    REAL(wp) :: clab, nlab, plab                        ! plant labile pools of carbon, nitrogen and phosphorus
    INTEGER  :: pftid                                   ! ID of the pft

    pftid = plant%pftid

    !> 1.0 IF: root and demand based N/P uptake   
    !! 
    !!
    IF(veg_config%flag_calc_soil_uptake) THEN ! case of root and demand based N/P uptake

      !> 1.1 calculate current "labile" C:N:P   
      !! 
      !!
      clab = plant%pool%leaf%carbon + plant%pool%fine_root%carbon + & 
             plant%pool%labile%carbon
      nlab = plant%pool%leaf%nitrogen + plant%pool%fine_root%nitrogen + & 
             plant%pool%labile%nitrogen
      plab = plant%pool%leaf%phosphorus + plant%pool%fine_root%phosphorus + & 
             plant%pool%labile%phosphorus

      !> 1.2 convert this into an instantaneous scalar of nutrient demand (above a target minimum C:N:P)
      !! 
      !! 
      IF(clab > eps8) THEN
         plant%f_n_demand = MAX(1._wp / veg_param%cn_leaf(pftid) - nlab / clab, 0.0_wp) 
      ELSE
         plant%f_n_demand = 1.0_wp
      ENDIF
      IF(nlab > eps8) THEN
         plant%f_p_demand = MAX(1._wp / veg_param%np_leaf(pftid) - plab / nlab, 0.0_wp)
      ELSE
         plant%f_p_demand = 1.0_wp
      ENDIF

      !> 1.3 average demand over a couple of days to get rid of diurnal C cycle:    
      !! ==> has been moved completely to tile%plant%moving_average -> PROCEDURE: update_soil
      !! 

      !> 1.4 calculate N, P uptake according to OCN 
      !! 
      !!
      plant%flux%uptake_n = veg_param%vmax_uptake_n(pftid) * &
              plant%moving_average%demand_uptake_n * & 
              soil%ammonium * ( veg_param%soil_n_aff1(pftid) + & 
                1._wp / ( veg_param%soil_n_aff2(pftid) + soil%ammonium )) * & 
              plant%pool%fine_root%carbon
      plant%flux%uptake_p = veg_param%vmax_uptake_p(pftid) * & 
              plant%moving_average%demand_uptake_p * & 
              soil%phosphate * ( veg_param%soil_p_aff1(pftid) + & 
                1._wp / ( veg_param%soil_p_aff2(pftid) + soil%phosphate )) * & 
              plant%pool%fine_root%carbon

    !> 2.0 ELSE: maintain prescribed C:N:P irrespective of nutrient uptake capacity == no nutrient limitation case 
    !! 
    !! @todo 2.0 is only a fallback option
    ELSE
      ! uptake required to maintain a presribed C:N ratio of the labile pool in µmol N / m2 / s
      plant%flux%uptake_n = MAX(plant%pool%labile%carbon / veg_param%cn_leaf(pftid) - & 
                                        plant%pool%labile%nitrogen, 0.0_wp) * 1.e6 / dt
      ! uptake required to maintain a presribed N:P ratio of the labile pool in µmol P / m2 / s
      plant%flux%uptake_p = MAX(plant%pool%labile%nitrogen / veg_param%np_leaf(pftid) + & 
                                           plant%flux%uptake_n * dt / 1.e6 - &
                                           plant%pool%labile%phosphorus, 0.0_wp) * 1.e6 / dt

    ENDIF     
  END SUBROUTINE update_vegsoil_bgc

  !-----------------------------------------------------------------------------------------------------
  !> update soil moisture and theta for the plant only model
  !! 
  !! two cases considered: a) dynamic water uptake calculation, b) prescribed values from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_physics(plant, soil)

    USE mo_base_memory,       ONLY: t_slm_soil
    USE mo_veg_memory,        ONLY: t_slm_plant
    USE mo_base_constants,    ONLY: dt
    USE mo_atmland_memory,    ONLY: atm2land
    USE mo_atmland_util,      ONLY: calc_spec_humidity_sat
    USE mo_land_constants,    ONLY: r_gas_dryair

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_plant), INTENT(inout) :: plant          !< plant variables per tile
    TYPE(t_slm_soil), INTENT(inout)  :: soil           !< soil variables per tile
    ! 0.2 Local
    INTEGER  :: t
    REAL(wp) :: q_sat
    REAL(wp) :: evaporation

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
