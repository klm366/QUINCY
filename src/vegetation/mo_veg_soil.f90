!!
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

  IMPLICIT NONE

  PUBLIC :: calc_nutrient_uptake

CONTAINS

  SUBROUTINE calc_nutrient_uptake(tile, veg_param)
    ! 
    ! routine to simulate N and P uptake from soil in the absence of a proper soil model
    ! - alternative use - supply N and P to maintain a certain target C:N:P of the labile pool
    USE mo_base_memory
    USE mo_base_util_functions, ONLY: moving_avrg
    USE mo_veg_parameters
    USE mo_veg_config_class,   ONLY: veg_config

    IMPLICIT NONE
    ! 0.1 Input
    TYPE(VegetationParameterClass),INTENT(in) :: veg_param     ! vegetation type parameters
    ! 0.2 Inout
    TYPE(TileClass), DIMENSION(ntile), INTENT(inout) :: tile
    ! 0.4 Local
    TYPE(AccumulationIntervalPlantClass) :: plant_acc_interval ! intervals over which variables are averaged
    REAL(wp) :: clab, nlab, plab ! plant labile pools of carbon, nitrogen and phosphorus
    REAL(wp) :: f_n_demand, f_p_demand ! relative plant demand for N and P 
    INTEGER :: t,pftid

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

  END SUBROUTINE calc_nutrient_uptake

END MODULE mo_veg_soil
