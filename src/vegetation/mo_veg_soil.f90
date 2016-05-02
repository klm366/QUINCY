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

  USE mo_kind, ONLY: wp
  USE mo_veg_config_class,   ONLY: veg_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_vegsoil_bgc,update_vegsoil_physics

CONTAINS

  SUBROUTINE update_vegsoil_bgc(tile, veg_param)
    ! 
    ! routine to simulate N and P uptake from soil in the absence of a proper soil model
    ! - alternative use - supply N and P to maintain a certain target C:N:P of the labile pool
    USE mo_base_memory
    USE mo_base_util,        ONLY: moving_avrg
    USE mo_veg_parameters

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_param),INTENT(in)                   :: veg_param     ! vegetation type parameters
    TYPE(t_slm_tile), DIMENSION(ntile), INTENT(inout)  :: tile
    ! 0.2 Local
    TYPE(AccumulationIntervalPlantClass) :: plant_acc_interval   ! intervals over which variables are averaged
    REAL(wp) :: clab, nlab, plab                                 ! plant labile pools of carbon, nitrogen and phosphorus
    REAL(wp) :: f_n_demand, f_p_demand                           ! relative plant demand for N and P 
    INTEGER :: t,pftid                                           ! iterate over tiles; ID of the pft

    DO t = 1,ntile
    
       pftid = tile(t)%plant%pftid

       IF(veg_config%flag_calc_soil_uptake) THEN ! case of root and demand based N/P uptake
          !
          ! 1. calculate current "labile" C:N:P   
          ! 
          clab = tile(t)%plant%pool%leaf%carbon + tile(t)%plant%pool%fine_root%carbon + & 
                 tile(t)%plant%pool%labile%carbon
          nlab = tile(t)%plant%pool%leaf%nitrogen + tile(t)%plant%pool%fine_root%nitrogen + & 
                 tile(t)%plant%pool%labile%nitrogen
          plab = tile(t)%plant%pool%leaf%phosphorus + tile(t)%plant%pool%fine_root%phosphorus + & 
                 tile(t)%plant%pool%labile%phosphorus
          !
          ! 2. convert this into an instantaneous scalar of nutrient demand (above a target minimum C:N:P)
          ! 
          IF(clab > eps8) THEN
             f_n_demand = MAX(1._wp / veg_param%cn_leaf(pftid) - nlab / clab, 0.0_wp) 
          ELSE
             f_n_demand = 1.0_wp
          ENDIF
          IF(nlab > eps8) THEN
             f_p_demand = MAX(1._wp / veg_param%np_leaf(pftid) - plab / nlab, 0.0_wp)
          ELSE
             f_p_demand = 1.0_wp
          ENDIF
          !
          ! 3. average demand over a couple of days to get rid of diurnal C cycle   
          ! 
          tile(t)%plant%accumulated%demand_uptake_n = moving_avrg( & 
                  tile(t)%plant%accumulated%demand_uptake_n, f_n_demand, & 
                  1.0_wp / (plant_acc_interval%demand_uptake * one_day/dt), &
                  tile(t)%plant%pool%fine_root%carbon < eps8 )
          tile(t)%plant%accumulated%demand_uptake_p = moving_avrg( & 
                  tile(t)%plant%accumulated%demand_uptake_p, f_p_demand, & 
                  1.0_wp / (plant_acc_interval%demand_uptake * one_day/dt), &
                  tile(t)%plant%pool%fine_root%carbon < eps8)
          ! 
          ! 4. calculate N, P uptake according to OCN  
          !
          tile(t)%plant%flux%uptake_n = veg_param%vmax_uptake_n(pftid) * &
                  tile(t)%plant%accumulated%demand_uptake_n * & 
                  tile(t)%soil%ammonium * ( veg_param%soil_n_aff1(pftid) + & 
                    1._wp / ( veg_param%soil_n_aff2(pftid) + tile(t)%soil%ammonium )) * & 
                  tile(t)%plant%pool%fine_root%carbon
          tile(t)%plant%flux%uptake_p = veg_param%vmax_uptake_p(pftid) * & 
                  tile(t)%plant%accumulated%demand_uptake_p * & 
                  tile(t)%soil%phosphate * ( veg_param%soil_p_aff1(pftid) + & 
                    1._wp / ( veg_param%soil_p_aff2(pftid) + tile(t)%soil%phosphate )) * & 
                  tile(t)%plant%pool%fine_root%carbon

        ELSE ! maintain prescribed C:N:P irrespective of nutrient uptake capacity == no nutrient limitation case
 
          ! uptake required to maintain a presribed C:N ratio of the labile pool in µmol N / m2 / s
          tile(t)%plant%flux%uptake_n = MAX(tile(t)%plant%pool%labile%carbon / veg_param%cn_leaf(pftid) - & 
                                            tile(t)%plant%pool%labile%nitrogen, 0.0_wp) * 1.e6 / dt
          ! uptake required to maintain a presribed N:P ratio of the labile pool in µmol P / m2 / s
          tile(t)%plant%flux%uptake_p = MAX(tile(t)%plant%pool%labile%nitrogen / veg_param%np_leaf(pftid) + & 
                                            tile(t)%plant%flux%uptake_n * dt / 1.e6 - &
                                            tile(t)%plant%pool%labile%phosphorus, 0.0_wp) * 1.e6 / dt

        ENDIF     
 
     ENDDO 

  END SUBROUTINE update_vegsoil_bgc

  PURE ELEMENTAL SUBROUTINE update_vegsoil_physics(tile)
    ! 
    ! Subroutine to update soil moisture and theta for the plant only model
    !
    USE mo_base_memory
    USE mo_atmland_memory
    USE mo_atmland_util, ONLY: calc_spec_humidity_sat
    USE mo_land_constants, ONLY: r_gas_dryair
    IMPLICIT NONE
    ! 0.2 Inout
    !TYPE(t_slm_tile), DIMENSION(ntile), INTENT(inout) :: tile
    TYPE(t_slm_tile), INTENT(inout) :: tile
    ! 0.4 local
    INTEGER :: t
    REAL(wp) :: q_sat
    REAL(wp) :: evaporation
  
    IF(veg_config%flag_calc_soil_moisture) THEN ! case of dynamic water uptake calculation 

       ! saturating specific humidity 

       q_sat = calc_spec_humidity_sat(atm2land%t_air,atm2land%press_srf)

       IF(tile%plant%gc > eps8)THEN
          ! T = rho * (qsat - qair ) * Gtot in mm/s 
          tile%plant%flux%transpiration = atm2land%press_srf / ( r_gas_dryair * atm2land%t_air ) * & 
                     (q_sat - atm2land%q_air) / (1._wp / tile%plant%gc + 1._wp / atm2land%ga )  
       ELSE 
          tile%plant%flux%transpiration = 0.0_wp 
       ENDIF
       ! E = rho * (qsat - qair) * theta * 1/ga_soil (taken as ga_soil = ga / 10)
       evaporation = atm2land%press_srf / ( r_gas_dryair * atm2land%t_air ) * & 
                     (q_sat - atm2land%q_air) * tile%soil%theta / 100.0_wp / ( 10._wp / atm2land%ga)

       ! update soil moisture dS = rain - transpiration - evaporation 
       tile%soil%moisture = tile%soil%moisture + (atm2land%rain + atm2land%snow & 
            - tile%plant%flux%transpiration - evaporation) * dt
       ! constrain to possible solutions
       tile%soil%moisture = MAX(MIN(tile%soil%moisture,veg_config%soil_awc_prescribe),0.0_wp)

       ! scale to relative water availability (%)
       tile%soil%theta = tile%soil%moisture / veg_config%soil_awc_prescribe * 100._wp

   ELSE 
       ! use prescribed value of theta from namelist
       tile%soil%theta = veg_config%soil_theta_prescribe 
 
       ! calculate corresponding soil moisture value (mm)     
       tile%soil%moisture = veg_config%soil_awc_prescribe * tile%soil%theta / 100.0_wp

   ENDIF 

  END SUBROUTINE update_vegsoil_physics

END MODULE mo_veg_soil
