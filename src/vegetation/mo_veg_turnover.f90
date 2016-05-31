!>
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
  USE mo_veg_memory

  PRIVATE

  PUBLIC :: calc_turnover

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Calculates litter fall resulting from tissue senesence, including nutrient resorption, 
  !! as well as turnover of sapwood ito heartwood and and fruits to seed_bed 
  !!
  !! @todo: mo_Veg_turnover: implement continuous turnover of foliar and root nutrients
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_turnover(plant,veg_param)

    USE mo_base_constants,     ONLY: one_day, one_year
    USE mo_shared_constants
    USE mo_veg_parameters

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_param), INTENT(in) :: veg_param 
    TYPE(t_slm_plant), INTENT(inout) :: plant
    ! 0.2 Local 
    TYPE(t_slm_veg_pool)                :: litterfall,growth
    REAL(wp)                        :: resorption_n,resorption_n15,resorption_p
    REAL(wp)                        :: fturn_leaf,fturn
    REAL(wp)                        :: nitrogen_resorption_fract, phosphorus_resorption_fract

    !> 1.0 initialise resporption and turnover rates
    !!
    !! 
    CALL growth%init(0.0_wp)
    CALL litterfall%init(0.0_wp)
    nitrogen_resorption_fract = update_nitrogen_resorption_fract(plant)
    phosphorus_resorption_fract = update_phosphorus_resorption_fract(plant)
    !> 2.0 foliage turnover 
    !! flat rate for evergreens, dual rates for deciduous
    !! includes retention of nutrients 
    SELECT CASE (veg_param%phenology_type(plant%pftid))
      CASE (ievergreen)
         fturn_leaf = veg_param%tau_leaf(plant%pftid) * dt/one_day/one_year 
      CASE (iraingreen, isummergreen)
         IF(plant%growing_season)THEN
            fturn_leaf = 0.0_wp              ! herbivory
         ELSE
            fturn_leaf = 0.2_wp * dt/one_day ! senescence
            ! if LAI dropy below minimum threshold, dump all leaves
            IF(plant%lai < min_lai) fturn_leaf = 1.0_wp
         ENDIF
      CASE (iperennial)
         IF(plant%growing_season)THEN
            fturn_leaf = veg_param%tau_leaf(plant%pftid) * dt/one_day/one_year
         ELSE
            fturn_leaf = 0.2_wp * dt/one_day
            IF(plant%lai < min_lai) fturn_leaf = 1.0_wp
       ENDIF
    END SELECT
    litterfall%leaf = fturn_leaf * plant%pool%leaf

    !> 2.1 Nutrient resportion from foliage turnover
    !! 
    resorption_n = nitrogen_resorption_fract * litterfall%leaf%nitrogen
    litterfall%leaf%nitrogen = litterfall%leaf%nitrogen - resorption_n 
    resorption_n15 = nitrogen_resorption_fract * litterfall%leaf%nitrogen15
    litterfall%leaf%nitrogen15 = litterfall%leaf%nitrogen15 - resorption_n15 
    resorption_p = phosphorus_resorption_fract * litterfall%leaf%phosphorus
    litterfall%leaf%phosphorus = litterfall%leaf%phosphorus - resorption_p

    !> 3.0 fine root turnover to litter
    !! 
    fturn = veg_param%tau_fine_root(plant%pftid) * dt/one_day/one_year
    litterfall%fine_root = fturn * plant%pool%fine_root

    !> 4.0 coarse root turnover to litter 
    !! 
    fturn = veg_param%tau_coarse_root(plant%pftid) * dt/one_day/one_year 
    litterfall%coarse_root = fturn * plant%pool%coarse_root

    !> 5.0 sapwood turnover 
    !> 5.1 litter production from halms (grasses), or branches (trees)
    !! 
    SELECT CASE (veg_param%growthform(plant%pftid))
      CASE (itree)
        fturn = veg_param%tau_branches(plant%pftid) * dt/one_day/one_year 
      CASE (igrass) ! for grases follows leaf phenology
        fturn = fturn_leaf
    ENDSELECT 
    litterfall%sap_wood = fturn * plant%pool%sap_wood

    !> 5.2 heart wood formation, including nutrient retention 
    !! 
    IF(veg_param%growthform(plant%pftid) == itree)THEN
        fturn = veg_param%tau_sap_wood(plant%pftid) * dt/one_day/one_year 
        growth%heart_wood = fturn * plant%pool%sap_wood
        growth%sap_wood = growth%sap_wood - growth%heart_wood 
        ! Nutrient retention 
        resorption_n = resorption_n + resorption_wood * growth%heart_wood%nitrogen 
        growth%heart_wood%nitrogen = (1._wp - resorption_wood) * growth%heart_wood%nitrogen
        resorption_n15 = resorption_n15 + resorption_wood * growth%heart_wood%nitrogen15 
        growth%heart_wood%nitrogen15 =  (1._wp - resorption_wood) * growth%heart_wood%nitrogen15
        resorption_p = resorption_p + phosphorus_resorption_fract * growth%heart_wood%phosphorus 
        growth%heart_wood%phosphorus = (1._wp - resorption_wood) * growth%heart_wood%phosphorus
    ENDIF 

    !> 6.0 fruit turnover to seed_bed
    !! 
    fturn = veg_param%tau_fruit(plant%pftid) * dt/one_day/one_year 
    growth%seed_bed = fturn * plant%pool%fruit
    growth%fruit = growth%fruit - growth%seed_bed

    !> 7.0 seed-bed turnover to litter 
    !! 
    fturn = veg_param%tau_seed(plant%pftid) * dt/one_day/one_year 
    litterfall%seed_bed = fturn * plant%pool%seed_bed

    !> 8.0 update litter, growth and resorption of the plant class
    !!
    plant%flux%litterfall     = plant%flux%litterfall     + litterfall
    plant%flux%growth         = plant%flux%growth         + growth
    plant%flux%resorption_n   = plant%flux%resorption_n   + resorption_n
    plant%flux%resorption_n15 = plant%flux%resorption_n15 + resorption_n15
    plant%flux%resorption_p   = plant%flux%resorption_p   + resorption_p

  END SUBROUTINE calc_turnover

  FUNCTION update_nitrogen_resorption_fract(plant) RESULT (resorption)
    !
    ! funtion to calculate fractional N resorption from leaves 
    !   should be made dynamically dependended on N stress later 
    !
    USE mo_veg_parameters

    IMPLICIT NONE
    TYPE(t_slm_plant) :: plant
    REAL(wp) :: resorption

    resorption = resorption_default

  END FUNCTION

  FUNCTION update_phosphorus_resorption_fract(plant) RESULT (resorption)
    !
    ! funtion to calculate fractional P resorption from leaves 
    !   should be made dynamically dependended on P stress later 
    !
    USE mo_veg_parameters

    IMPLICIT NONE
    TYPE(t_slm_plant) :: plant
    REAL(wp) :: resorption 

    resorption = resorption_default

  END FUNCTION

END MODULE mo_veg_turnover
