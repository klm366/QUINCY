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
  USE mo_jsb_lctlib_class,     ONLY: t_lctlib_element
  USE mo_veg_memory
  USE mo_veg_constants


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_turnover

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Calculates litter fall resulting from tissue senesence, including nutrient recycling, 
  !! as well as turnover of sapwood to heartwood and fruits to seed_bed 
  !!
  !! @todo: mo_veg_turnover:calc_turnover -> herbivory on leaves as turnover term? 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_turnover(tile, lctlib)

    USE mo_base_memory,        ONLY: t_slm_tile
    USE mo_base_constants,     ONLY: one_day, one_year
    USE mo_shared_constants

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile      !< one tile
    TYPE(t_lctlib_element),INTENT(in)       :: lctlib    !< land-cover-type library - parameters for one pft
    ! 0.2 Local 
    TYPE(t_slm_plant),POINTER         :: plant     ! plant variables per tile
    TYPE(t_slm_veg_pool)              :: litterfall,growth
    REAL(wp)                          :: recycling_leaf_n,recycling_leaf_n15,recycling_leaf_p
    REAL(wp)                          :: recycling_sap_wood_n,recycling_sap_wood_n15,recycling_sap_wood_p
    REAL(wp)                          :: fturn_leaf,fturn
    REAL(wp)                          :: nitrogen_recycling_fract, phosphorus_recycling_fract

    ! bi-directional pointer
    plant => tile%plant


    !> 1.0 initialise output variables
    !!
    CALL growth%init(0.0_wp)
    CALL litterfall%init(0.0_wp)
    recycling_leaf_n = 0.0_wp
    recycling_leaf_p = 0.0_wp
    recycling_leaf_n15 = 0.0_wp
    recycling_sap_wood_n = 0.0_wp
    recycling_sap_wood_p = 0.0_wp
    recycling_sap_wood_n15 = 0.0_wp

    !> 2.0 foliage turnover, including recycling of nutrients 
    !!
    !! flat rate for evergreens, dual rates for deciduous (senescence after the growing season)
    SELECT CASE (lctlib%phenology_type)
      CASE (ievergreen)
         fturn_leaf = lctlib%tau_leaf * dt/one_day/one_year 
      CASE (iraingreen, isummergreen)
         IF(plant%growing_season)THEN
            fturn_leaf = 0.0_wp              ! herbivory
         ELSE
            fturn_leaf = 0.2_wp * dt/one_day ! senescence
            ! if LAI drops below minimum threshold, dump all leaves
            IF(plant%lai < min_lai) fturn_leaf = 1.0_wp
         ENDIF
      CASE (iperennial)
         IF(plant%growing_season)THEN
            fturn_leaf = lctlib%tau_leaf * dt/one_day/one_year
         ELSE
            ! if LAI drops below minimum threshold, dump all leaves
            fturn_leaf = 0.2_wp * dt/one_day
            IF(plant%lai < min_lai) fturn_leaf = 1.0_wp
       ENDIF
    END SELECT
    litterfall%leaf = fturn_leaf * plant%pool%leaf

    !> 2.1 Nutrient recycling from foliage turnover
    !! 
    nitrogen_recycling_fract = update_nitrogen_recycling_fract(plant)
    phosphorus_recycling_fract = update_phosphorus_recycling_fract(plant)
    recycling_leaf_n = nitrogen_recycling_fract * litterfall%leaf%nitrogen
    litterfall%leaf%nitrogen = litterfall%leaf%nitrogen - recycling_leaf_n 
    recycling_leaf_n15 = nitrogen_recycling_fract * litterfall%leaf%nitrogen15
    litterfall%leaf%nitrogen15 = litterfall%leaf%nitrogen15 - recycling_leaf_n15 
    recycling_leaf_p = phosphorus_recycling_fract * litterfall%leaf%phosphorus
    litterfall%leaf%phosphorus = litterfall%leaf%phosphorus - recycling_leaf_p

    !> 3.0 fine root turnover to litter
    !! 
    fturn = lctlib%tau_fine_root * dt/one_day/one_year
    litterfall%fine_root = fturn * plant%pool%fine_root

    !> 4.0 coarse root turnover to litter 
    !! 
    fturn = lctlib%tau_coarse_root * dt/one_day/one_year 
    litterfall%coarse_root = fturn * plant%pool%coarse_root

    !> 5.0 sapwood turnover 
    !! 
    !> 5.1 litter production from halms (grasses), or branches (trees)
    SELECT CASE (lctlib%growthform)
      CASE (itree)
        fturn = lctlib%tau_branches * dt/one_day/one_year 
      CASE (igrass) ! for grases follows leaf phenology
        fturn = fturn_leaf
    ENDSELECT 
    litterfall%sap_wood = fturn * plant%pool%sap_wood

    !> 5.2 heart wood formation, including nutrient recycling
    !! 
    IF(lctlib%growthform == itree)THEN
        fturn = lctlib%tau_sap_wood * dt/one_day/one_year 
        growth%heart_wood = fturn * plant%pool%sap_wood
        growth%sap_wood = growth%sap_wood - growth%heart_wood 
        ! Nutrient recycling
        nitrogen_recycling_fract = update_nitrogen_recycling_fract(plant)
        phosphorus_recycling_fract = update_phosphorus_recycling_fract(plant)
        recycling_sap_wood_n = nitrogen_recycling_fract * growth%heart_wood%nitrogen 
        growth%heart_wood%nitrogen = (1._wp - nitrogen_recycling_fract) * growth%heart_wood%nitrogen
        recycling_sap_wood_n15 = nitrogen_recycling_fract * growth%heart_wood%nitrogen15 
        growth%heart_wood%nitrogen15 =  (1._wp - nitrogen_recycling_fract) * growth%heart_wood%nitrogen15
        recycling_sap_wood_p = phosphorus_recycling_fract * growth%heart_wood%phosphorus 
        growth%heart_wood%phosphorus = (1._wp - phosphorus_recycling_fract) * growth%heart_wood%phosphorus
    ENDIF 

    !> 6.0 fruit turnover to seed_bed
    !! 
    fturn = lctlib%tau_fruit * dt/one_day/one_year 
    growth%seed_bed = fturn * plant%pool%fruit
    growth%fruit = growth%fruit - growth%seed_bed

    !> 7.0 seed-bed turnover to litter 
    !! 
    fturn = lctlib%tau_seed_litter * dt/one_day/one_year 
    litterfall%seed_bed = fturn * plant%pool%seed_bed

    !> 8.0 update litter, growth and recycling of the plant class
    !!
    plant%litterfall             = plant%litterfall             + litterfall
    plant%growth                 = plant%growth                 + growth
    plant%recycling_leaf_n       = plant%recycling_leaf_n       + recycling_leaf_n
    plant%recycling_leaf_n15     = plant%recycling_leaf_n15     + recycling_leaf_n15
    plant%recycling_leaf_p       = plant%recycling_leaf_p       + recycling_leaf_p
    plant%recycling_sap_wood_n   = plant%recycling_sap_wood_n   + recycling_sap_wood_n
    plant%recycling_sap_wood_n15 = plant%recycling_sap_wood_n15 + recycling_sap_wood_n15
    plant%recycling_sap_wood_p   = plant%recycling_sap_wood_p   + recycling_sap_wood_p

  END SUBROUTINE calc_turnover

  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional N recycling from leaves 
  !!
  !! should be made dynamically dependended on N stress later 
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION update_nitrogen_recycling_fract(plant) RESULT (recycling)

    IMPLICIT NONE

    TYPE(t_slm_plant) :: plant
    REAL(wp) :: recycling

    recycling = resorption_default

  END FUNCTION

  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional P recycling from leaves
  !!
  !! should be made dynamically dependended on P stress later
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION update_phosphorus_recycling_fract(plant) RESULT (recycling)

    IMPLICIT NONE

    TYPE(t_slm_plant) :: plant
    REAL(wp) :: recycling 

    recycling = resorption_default

  END FUNCTION

END MODULE mo_veg_turnover
