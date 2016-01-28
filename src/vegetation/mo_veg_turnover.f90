!!
!! @brief Contains the routines for plant turnover related processes 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_turnover

  USE mo_kind
  USE mo_base_memory
  USE mo_veg_memory

  PRIVATE

  REAL(wp) :: fturn_leaf,fturn_fine_root,fturn_coarse_root,fturn_sap_wood, &
                fturn_heart_wood,fturn_seed_bed 

  PUBLIC :: calc_turnover

CONTAINS

  SUBROUTINE calc_turnover(plant,veg_param)
    !
    ! Routine to calculate litter fall resulting from tissue senesence, 
    ! as well as turnover of sapwood and fruit pools,
    ! and turnover of foliar and root nutrients (not yet implemented) 
    ! 
    USE mo_land_constants
    USE mo_veg_parameters

    IMPLICIT NONE
    ! 0.1 Input
    TYPE(VegetationParameterClass), INTENT(in) :: veg_param 
    ! 0.2 In- and Output
    TYPE(PlantClass), INTENT(inout) :: plant
    ! 0.4 Local 
    REAL(wp)                        :: nitrogen_resorption, phosphorus_resorption

    !
    ! 1.0 initialise resporption and turnover rates
    !
    nitrogen_resorption = update_nitrogen_resorption(plant)
    phosphorus_resorption = update_phosphorus_resorption(plant)
    CALL update_turnover_fraction(plant)

    !
    ! 2.1 C Litter fall
    !
    plant%flux%litterfall%leaf%carbon = fturn_leaf * plant%pool%leaf%carbon 
    plant%flux%litterfall%fine_root%carbon = fturn_fine_root * plant%pool%fine_root%carbon 
    plant%flux%litterfall%coarse_root%carbon = fturn_coarse_root * plant%pool%coarse_root%carbon 
    plant%flux%litterfall%sap_wood%carbon = fturn_sap_wood * plant%pool%sap_wood%carbon 
    plant%flux%litterfall%heart_wood%carbon = fturn_heart_wood * plant%pool%heart_wood%carbon 
    plant%flux%litterfall%seed_bed%carbon = fturn_seed_bed * plant%pool%seed_bed%carbon 

    !
    ! 2.2 N Litter fall: a fraction of foliar N becomes resorpted into the labile pool
    ! not entirely happy with this solution of resportion. Ideally the resportion should be 
    ! calculated as a separate function and costly in terms of C
    !
    plant%flux%litterfall%leaf%nitrogen = fturn_leaf * (1.0_wp -nitrogen_resorption) * plant%pool%leaf%nitrogen 
    plant%flux%growth%labile%nitrogen = plant%flux%growth%labile%nitrogen + & 
                                        fturn_leaf * nitrogen_resorption * plant%pool%leaf%nitrogen 
    plant%flux%growth%leaf%nitrogen = plant%flux%growth%leaf%nitrogen - & 
                                        fturn_leaf * nitrogen_resorption * plant%pool%leaf%nitrogen 
    plant%flux%litterfall%fine_root%nitrogen = fturn_fine_root * plant%pool%fine_root%nitrogen 
    plant%flux%litterfall%coarse_root%nitrogen = fturn_coarse_root * plant%pool%coarse_root%nitrogen 
    plant%flux%litterfall%sap_wood%nitrogen = fturn_sap_wood * plant%pool%sap_wood%nitrogen 
    plant%flux%litterfall%heart_wood%nitrogen = fturn_heart_wood * plant%pool%heart_wood%nitrogen 
    plant%flux%litterfall%seed_bed%nitrogen = fturn_seed_bed * plant%pool%seed_bed%nitrogen 

    !
    ! 2.3 P Litter fall: a fraction of foliar P becomes resorpted into the labile pool
    ! not entirely happy with this solution of resportion. Ideally the resportion should be 
    ! calculated as a separate function and costly in terms of C
    !
    plant%flux%litterfall%leaf%phosphorus = fturn_leaf * (1.0_wp-phosphorus_resorption) * plant%pool%leaf%phosphorus 
    plant%flux%growth%labile%phosphorus = plant%flux%growth%labile%phosphorus + & 
                                        fturn_leaf * phosphorus_resorption * plant%pool%leaf%phosphorus 
    plant%flux%growth%leaf%phosphorus = plant%flux%growth%leaf%phosphorus - & 
                                        fturn_leaf * phosphorus_resorption * plant%pool%leaf%phosphorus 
    plant%flux%litterfall%fine_root%phosphorus = fturn_fine_root * plant%pool%fine_root%phosphorus 
    plant%flux%litterfall%coarse_root%phosphorus = fturn_coarse_root * plant%pool%coarse_root%phosphorus 
    plant%flux%litterfall%sap_wood%phosphorus = fturn_sap_wood * plant%pool%sap_wood%phosphorus 
    plant%flux%litterfall%heart_wood%phosphorus = fturn_heart_wood * plant%pool%heart_wood%phosphorus 
    plant%flux%litterfall%seed_bed%phosphorus = fturn_seed_bed * plant%pool%seed_bed%phosphorus 

    !
    ! 2.4 C13 Litter fall
    !
    plant%flux%litterfall%leaf%carbon13 = fturn_leaf * plant%pool%leaf%carbon13 
    plant%flux%litterfall%fine_root%carbon13 = fturn_fine_root * plant%pool%fine_root%carbon13 
    plant%flux%litterfall%coarse_root%carbon13 = fturn_coarse_root * plant%pool%coarse_root%carbon13 
    plant%flux%litterfall%sap_wood%carbon13 = fturn_sap_wood * plant%pool%sap_wood%carbon13 
    plant%flux%litterfall%heart_wood%carbon13 = fturn_heart_wood * plant%pool%heart_wood%carbon13 
    plant%flux%litterfall%seed_bed%carbon13 = fturn_seed_bed * plant%pool%seed_bed%carbon13 

    !
    ! 2.5 C14 Litter fall
    !
    plant%flux%litterfall%leaf%carbon14 = fturn_leaf * plant%pool%leaf%carbon14 
    plant%flux%litterfall%fine_root%carbon14 = fturn_fine_root * plant%pool%fine_root%carbon14 
    plant%flux%litterfall%coarse_root%carbon14 = fturn_coarse_root * plant%pool%coarse_root%carbon14 
    plant%flux%litterfall%sap_wood%carbon14 = fturn_sap_wood * plant%pool%sap_wood%carbon14 
    plant%flux%litterfall%heart_wood%carbon14 = fturn_heart_wood * plant%pool%heart_wood%carbon14 
    plant%flux%litterfall%seed_bed%carbon14 = fturn_seed_bed * plant%pool%seed_bed%carbon14 

    !
    ! 2.6 N15 Litter fall
    !
    plant%flux%litterfall%leaf%nitrogen15 = fturn_leaf * (1.0_wp-nitrogen_resorption) * plant%pool%leaf%nitrogen15 
    plant%flux%growth%labile%nitrogen15 = plant%flux%growth%labile%nitrogen15 + & 
                                        fturn_leaf * nitrogen_resorption * plant%pool%leaf%nitrogen15 
    plant%flux%growth%leaf%nitrogen15 = plant%flux%growth%leaf%nitrogen15 - & 
                                        fturn_leaf * nitrogen_resorption * plant%pool%leaf%nitrogen15
    plant%flux%litterfall%fine_root%nitrogen15 = fturn_fine_root * plant%pool%fine_root%nitrogen15 
    plant%flux%litterfall%coarse_root%nitrogen15 = fturn_coarse_root * plant%pool%coarse_root%nitrogen15 
    plant%flux%litterfall%sap_wood%nitrogen15 = fturn_sap_wood * plant%pool%sap_wood%nitrogen15 
    plant%flux%litterfall%heart_wood%nitrogen15 = fturn_heart_wood * plant%pool%heart_wood%nitrogen15 
    plant%flux%litterfall%seed_bed%nitrogen15 = fturn_seed_bed * plant%pool%seed_bed%nitrogen15 

    !
    ! calculate transformation from sapwood into heart wood 
    ! 
    CALL calc_pool_turnover(plant%pool%sap_wood,veg_param%tau_sapwood(plant%pftid),resorption_wood, & 
                            plant%flux%growth%sap_wood,plant%flux%growth%heart_wood,plant%flux%growth%labile)
    !
    ! calculate transformation from fruit into seed_bed 
    ! 
    CALL calc_pool_turnover(plant%pool%fruit,veg_param%tau_fruit_pool(plant%pftid),zero, & 
                            plant%flux%growth%fruit,plant%flux%growth%seed_bed,plant%flux%growth%labile)

  END SUBROUTINE calc_turnover

  SUBROUTINE update_turnover_fraction(plant)
    !
    ! routine to calculate fractional loss of biomass from a plant compartment
    ! set to constant for now
    ! should be based on PFT, phenology and senesence
    !
    USE mo_base_constants, ONLY: dt, one_day, one_year

    IMPLICIT NONE

    TYPE(PlantClass), INTENT(inout) :: plant

    fturn_leaf = 0.3_wp * dt/one_day/one_year
    fturn_fine_root = 0.7_wp * dt/one_day/one_year
    fturn_coarse_root = 0.1_wp * dt/one_day/one_year
    fturn_sap_wood = 0.01_wp * dt/one_day/one_year
    fturn_heart_wood = 0.01_wp * dt/one_day/one_year
    fturn_seed_bed = 1.0_wp * dt/one_day/one_year

  END SUBROUTINE update_turnover_fraction

  SUBROUTINE calc_pool_turnover(source_pool, tau_pool, resorption, source_flux, sink_flux, storage_flux) 
    !
    ! calculates mass transfer from one pool to another due to senesences
    !   used for sapwood -> heartwood
    !            fruit -> seed_bed 
    !   this is the only process forming these pools, added to the growth term of the respective 
    !   sink pools
    ! 
    USE mo_base_constants, ONLY: dt, one_day, one_year
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: tau_pool
    REAL(wp), INTENT(in) :: resorption

    TYPE(Element),INTENT(inout) :: source_flux,sink_flux,storage_flux,source_pool

    source_flux%carbon = source_flux%carbon - tau_pool * dt/one_day/one_year * source_pool%carbon
    sink_flux%carbon = tau_pool * dt/one_day/one_year * source_pool%carbon

    source_flux%nitrogen = source_flux%nitrogen - tau_pool * dt/one_day/one_year * source_pool%nitrogen
    sink_flux%nitrogen = ( 1.0_wp - resorption ) * tau_pool * dt/one_day/one_year * source_pool%nitrogen
    storage_flux%nitrogen = storage_flux%nitrogen + resorption * tau_pool * dt/one_day/one_year * source_pool%nitrogen

    source_flux%phosphorus = source_flux%phosphorus - tau_pool * dt/one_day/one_year * source_pool%phosphorus
    sink_flux%phosphorus = ( 1.0_wp - resorption )  * tau_pool * dt/one_day/one_year * source_pool%phosphorus
    storage_flux%phosphorus = storage_flux%phosphorus + resorption * tau_pool * dt/one_day/one_year * source_pool%phosphorus

    source_flux%carbon13 = source_flux%carbon13 - tau_pool * dt/one_day/one_year * source_pool%carbon13
    sink_flux%carbon13 = tau_pool * dt/one_day/one_year * source_pool%carbon13

    source_flux%carbon14 = source_flux%carbon14 - tau_pool * dt/one_day/one_year * source_pool%carbon14
    sink_flux%carbon14 = tau_pool * dt/one_day/one_year * source_pool%carbon14

    source_flux%nitrogen15 = source_flux%nitrogen15 - tau_pool * dt/one_day/one_year * source_pool%nitrogen15
    sink_flux%nitrogen15 = ( 1.0_wp - resorption ) * tau_pool * dt/one_day/one_year * source_pool%nitrogen15
    storage_flux%nitrogen15 = storage_flux%nitrogen15 + resorption * tau_pool * dt/one_day/one_year * source_pool%nitrogen15

  END SUBROUTINE calc_pool_turnover

  FUNCTION update_nitrogen_resorption(plant) RESULT (resorption)
    !
    ! funtion to calculate fractional N resorption from leaves 
    !   should be made dynamically dependended on N stress later 
    !
    USE mo_veg_parameters

    IMPLICIT NONE
    TYPE(PlantClass) :: plant
    REAL(wp) :: resorption

    resorption = resorption_default

  END FUNCTION

  FUNCTION update_phosphorus_resorption(plant) RESULT (resorption)
    !
    ! funtion to calculate fractional P resorption from leaves 
    !   should be made dynamically dependended on P stress later 
    !
    USE mo_veg_parameters

    IMPLICIT NONE
    TYPE(PlantClass) :: plant
    REAL(wp) :: resorption 

    resorption = resorption_default

  END FUNCTION

END MODULE mo_veg_turnover
