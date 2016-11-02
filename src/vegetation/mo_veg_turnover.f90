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

#include "dsl_jsbach.inc"

  USE mo_kind

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
  SUBROUTINE calc_turnover(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_base_constants,        ONLY: dt, one_day, one_year
    USE mo_shared_constants
    USE mo_veg_constants

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout) :: tile      !< one tile
    ! 0.2 Local 
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp)                          :: recycling_leaf_n, &
                                         recycling_leaf_n15, &
                                         recycling_leaf_p
    REAL(wp)                          :: recycling_sap_wood_n, &
                                         recycling_sap_wood_n15, &
                                         recycling_sap_wood_p
    REAL(wp)                          :: fturn_leaf, &
                                         fturn
    REAL(wp)                          :: nitrogen_recycling_fract, &
                                         phosphorus_recycling_fract
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(VEG_)%p) 
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    veg    => veg__mem 
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

    !> 1.0 initialise output variables
    !!
    CALL veg%growth_t%init(0.0_wp) 
    CALL veg%litterfall_t%init(0.0_wp)
    recycling_leaf_n        = 0.0_wp
    recycling_leaf_n15      = 0.0_wp
    recycling_leaf_p        = 0.0_wp
    recycling_sap_wood_n    = 0.0_wp
    recycling_sap_wood_n15  = 0.0_wp
    recycling_sap_wood_p    = 0.0_wp

    !> 2.0 foliage turnover, including recycling of nutrients 
    !!
    !! flat rate for evergreens, dual rates for deciduous (senescence after the growing season)
    SELECT CASE (lctlib%phenology_type)
      CASE (ievergreen)
         fturn_leaf = lctlib%tau_leaf * dt/one_day/one_year 
      CASE (iraingreen, isummergreen)
         IF(veg%growing_season%ptr)THEN
            fturn_leaf = 0.0_wp              ! herbivory
         ELSE
            fturn_leaf = 0.2_wp * dt/one_day ! senescence
            ! if LAI drops below minimum threshold, dump all leaves
            IF(veg%lai%ptr < min_lai) fturn_leaf = 1.0_wp
         ENDIF
      CASE (iperennial)
         IF(veg%growing_season%ptr)THEN
            fturn_leaf = lctlib%tau_leaf * dt/one_day/one_year
         ELSE
            ! if LAI drops below minimum threshold, dump all leaves
            fturn_leaf = 0.2_wp * dt/one_day
            IF(veg%lai%ptr < min_lai) fturn_leaf = 1.0_wp
       ENDIF
    END SELECT
    veg%litterfall_t%leaf = veg%pool%leaf * fturn_leaf

    !> 2.1 Nutrient recycling from foliage turnover
    !! 
    nitrogen_recycling_fract    = update_nitrogen_recycling_fract()
    phosphorus_recycling_fract  = update_phosphorus_recycling_fract()

    recycling_leaf_n                = nitrogen_recycling_fract * veg%litterfall_t%leaf%nitrogen%ptr
    veg%litterfall_t%leaf%nitrogen    = veg%litterfall_t%leaf%nitrogen - recycling_leaf_n 
    recycling_leaf_n15              = nitrogen_recycling_fract * veg%litterfall_t%leaf%nitrogen15%ptr
    veg%litterfall_t%leaf%nitrogen15  = veg%litterfall_t%leaf%nitrogen15 - recycling_leaf_n15 
    recycling_leaf_p                = phosphorus_recycling_fract * veg%litterfall_t%leaf%phosphorus%ptr
    veg%litterfall_t%leaf%phosphorus  = veg%litterfall_t%leaf%phosphorus - recycling_leaf_p

    !> 3.0 fine root turnover to litter
    !! 
    fturn = lctlib%tau_fine_root * dt/one_day/one_year
    veg%litterfall_t%fine_root = veg%pool%fine_root * fturn

    !> 4.0 coarse root turnover to litter 
    !! 
    fturn = lctlib%tau_coarse_root * dt/one_day/one_year 
    veg%litterfall_t%coarse_root = veg%pool%coarse_root * fturn

    !> 5.0 sapwood turnover 
    !! 
    !> 5.1 litter production from halms (grasses), or branches (trees)
    SELECT CASE (lctlib%growthform)
      CASE (itree)
        fturn = lctlib%tau_branches * dt/one_day/one_year 
      CASE (igrass) ! for grases follows leaf phenology
        fturn = fturn_leaf
    ENDSELECT 
    veg%litterfall_t%sap_wood = veg%pool%sap_wood * fturn

    !> 5.2 heart wood formation, including nutrient recycling
    !! 
    IF(lctlib%growthform == itree)THEN
        fturn = lctlib%tau_sap_wood * dt/one_day/one_year 
        veg%growth_t%heart_wood = veg%pool%sap_wood     * fturn
        veg%growth_t%sap_wood   = veg%growth_t%sap_wood - veg%growth_t%heart_wood 

        ! Nutrient recycling
        nitrogen_recycling_fract    = update_nitrogen_recycling_fract()
        phosphorus_recycling_fract  = update_phosphorus_recycling_fract()
        recycling_sap_wood_n              = veg%growth_t%heart_wood%nitrogen%ptr   * nitrogen_recycling_fract
        veg%growth_t%heart_wood%nitrogen    = veg%growth_t%heart_wood%nitrogen       * (1._wp - nitrogen_recycling_fract)
        recycling_sap_wood_n15            = veg%growth_t%heart_wood%nitrogen15%ptr * nitrogen_recycling_fract
        veg%growth_t%heart_wood%nitrogen15  = veg%growth_t%heart_wood%nitrogen15     * (1._wp - nitrogen_recycling_fract)
        recycling_sap_wood_p              = veg%growth_t%heart_wood%phosphorus%ptr * phosphorus_recycling_fract
        veg%growth_t%heart_wood%phosphorus  = veg%growth_t%heart_wood%phosphorus     * (1._wp - phosphorus_recycling_fract)
    ENDIF 

    !> 6.0 fruit turnover to seed_bed
    !! 
    fturn = lctlib%tau_fruit * dt/one_day/one_year 
    veg%growth_t%seed_bed = veg%pool%fruit     * fturn
    veg%growth_t%fruit    = veg%growth_t%fruit - veg%growth_t%seed_bed

    !> 7.0 seed-bed turnover to litter 
    !! 
    fturn = lctlib%tau_seed_litter * dt/one_day/one_year 
    veg%litterfall_t%seed_bed = veg%pool%seed_bed * fturn

    !> 8.0 update litter, growth and recycling of the plant class
    !!
    veg%litterfall             = veg%litterfall + veg%litterfall_t
    veg%growth                 = veg%growth     + veg%growth_t

    veg%recycling_leaf_n       = veg%recycling_leaf_n       + recycling_leaf_n
    veg%recycling_leaf_n15     = veg%recycling_leaf_n15     + recycling_leaf_n15
    veg%recycling_leaf_p       = veg%recycling_leaf_p       + recycling_leaf_p
    veg%recycling_sap_wood_n   = veg%recycling_sap_wood_n   + recycling_sap_wood_n
    veg%recycling_sap_wood_n15 = veg%recycling_sap_wood_n15 + recycling_sap_wood_n15
    veg%recycling_sap_wood_p   = veg%recycling_sap_wood_p   + recycling_sap_wood_p

  END SUBROUTINE calc_turnover

  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional N recycling from leaves 
  !!
  !! should be made dynamically dependended on N stress later 
  !!
  !! @todo veg_turnover:update_nitrogen_recycling_fract implement functionality other than "set to default"
  !-----------------------------------------------------------------------------------------------------
  FUNCTION update_nitrogen_recycling_fract() RESULT (recycling)

    USE mo_veg_constants,     ONLY: resorption_default

    IMPLICIT NONE
    ! 0.1 InOut

    ! 0.2 Local
    REAL(wp) :: recycling


    recycling = resorption_default

  END FUNCTION

  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional P recycling from leaves
  !!
  !! should be made dynamically dependended on P stress later
  !!
  !! @todo veg_turnover:update_phosphorus_recycling_fract implement functionality other than "set to default"
  !-----------------------------------------------------------------------------------------------------
  FUNCTION update_phosphorus_recycling_fract() RESULT (recycling)

    USE mo_veg_constants,     ONLY: resorption_default

    IMPLICIT NONE
    ! 0.1 InOut

    ! 0.2 Local
    REAL(wp) :: recycling 


    recycling = resorption_default

  END FUNCTION

END MODULE mo_veg_turnover
