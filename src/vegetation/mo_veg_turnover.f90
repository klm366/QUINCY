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
  PUBLIC :: update_veg_turnover

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> Calculates litter fall resulting from tissue senesence, including nutrient recycling, 
  !! as well as turnover of sapwood to heartwood and fruits to seed_bed 
  !!
  !! @todo: mo_veg_turnover:update_veg_turnover -> herbivory on leaves as turnover term? 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_turnover(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout) :: tile      !< one tile
    ! 0.2 Local 
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)

    veg    => veg__mem 
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

    !> Sub-Task
    !! 
    CALL calc_turnover(lctlib, &
                        veg%pool, veg%growth, veg%litterfall, &
                        veg%growth_t, veg%litterfall_t, &
                        veg%growing_season%ptr, veg%lai%ptr, &
                        veg%recycling_leaf_n%ptr, veg%recycling_leaf_n15%ptr, veg%recycling_leaf_p%ptr, &
                        veg%recycling_sap_wood_n%ptr, veg%recycling_sap_wood_n15%ptr, veg%recycling_sap_wood_p%ptr)


  END SUBROUTINE update_veg_turnover


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_turnover
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover: recycling and retranslocation 
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_turnover(lctlib                     , &
!  ELEMENTAL SUBROUTINE calc_turnover(lctlib                     , &
                                      veg_pool                  , &
                                      veg_growth                , &
                                      veg_litterfall            , &
                                      veg_growth_t              , &
                                      veg_litterfall_t          , &
                                      growing_season            , &
                                      lai                       , &
                                      recycling_leaf_n          , &
                                      recycling_leaf_n15        , &
                                      recycling_leaf_p          , &
                                      recycling_sap_wood_n      , &
                                      recycling_sap_wood_n15    , &
                                      recycling_sap_wood_p      )

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_memory_class,      ONLY: t_slm_struct_veg_pool
    USE mo_base_constants,        ONLY: dt, one_day, one_year
    USE mo_shared_constants
    USE mo_veg_constants

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_lctlib_element),       INTENT(in)    :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_slm_struct_veg_pool),  INTENT(inout) :: veg_pool         , & !< 
                                                   veg_growth       , & !< 
                                                   veg_litterfall   , & !< 
                                                   veg_growth_t     , & !< 
                                                   veg_litterfall_t     !< 
    LOGICAL,                      INTENT(in)    :: growing_season       !< 
    REAL(wp),                     INTENT(in)    :: lai                  !< 
    REAL(wp),                     INTENT(inout) :: recycling_leaf_n          , & !< 
                                                   recycling_leaf_n15        , & !< 
                                                   recycling_leaf_p          , & !< 
                                                   recycling_sap_wood_n      , & !< 
                                                   recycling_sap_wood_n15    , & !< 
                                                   recycling_sap_wood_p          !< 
    ! 0.2 Local
    REAL(wp)                          :: retranslocation_leaf_n       , & !< leaf
                                         retranslocation_leaf_n15     , &
                                         retranslocation_leaf_p
    REAL(wp)                          :: retranslocation_sap_wood_n   , & !< sap_wood
                                         retranslocation_sap_wood_n15 , &
                                         retranslocation_sap_wood_p
    REAL(wp)                          :: fturn_leaf                   , & !< fturn
                                         fturn
    REAL(wp)                          :: nitrogen_recycling_fract     , & !< recycling fraction
                                         phosphorus_recycling_fract


    !> 1.0 initialise variables
    !!
    CALL veg_growth_t%init(0.0_wp) 
    CALL veg_litterfall_t%init(0.0_wp)
    retranslocation_leaf_n        = 0.0_wp
    retranslocation_leaf_n15      = 0.0_wp
    retranslocation_leaf_p        = 0.0_wp
    retranslocation_sap_wood_n    = 0.0_wp
    retranslocation_sap_wood_n15  = 0.0_wp
    retranslocation_sap_wood_p    = 0.0_wp

    !> 2.0 foliage turnover, including recycling of nutrients 
    !!
    !! flat rate for evergreens, dual rates for deciduous (senescence after the growing season)
    SELECT CASE (lctlib%phenology_type)
      CASE (ievergreen)
         fturn_leaf = lctlib%tau_leaf * dt/one_day/one_year 
      CASE (iraingreen, isummergreen)
         IF(growing_season)THEN
            fturn_leaf = 0.0_wp              ! herbivory
         ELSE
            fturn_leaf = 0.2_wp * dt/one_day ! senescence
            ! if LAI drops below minimum threshold, dump all leaves
            IF(lai < min_lai) fturn_leaf = 1.0_wp
         ENDIF
      CASE (iperennial)
         IF(growing_season)THEN
            fturn_leaf = lctlib%tau_leaf * dt/one_day/one_year
         ELSE
            ! if LAI drops below minimum threshold, dump all leaves
            fturn_leaf = 0.2_wp * dt/one_day
            IF(lai < min_lai) fturn_leaf = 1.0_wp
       ENDIF
    END SELECT
    veg_litterfall_t%leaf = veg_pool%leaf * fturn_leaf

    !> 2.1 Nutrient recycling from foliage turnover
    !! 
    nitrogen_recycling_fract    = calc_nitrogen_recycling_fract()
    phosphorus_recycling_fract  = calc_phosphorus_recycling_fract()

    retranslocation_leaf_n            = nitrogen_recycling_fract          * veg_litterfall_t%leaf%nitrogen%ptr
    veg_litterfall_t%leaf%nitrogen    = veg_litterfall_t%leaf%nitrogen    - retranslocation_leaf_n 
    retranslocation_leaf_n15          = nitrogen_recycling_fract          * veg_litterfall_t%leaf%nitrogen15%ptr
    veg_litterfall_t%leaf%nitrogen15  = veg_litterfall_t%leaf%nitrogen15  - retranslocation_leaf_n15 
    retranslocation_leaf_p            = phosphorus_recycling_fract        * veg_litterfall_t%leaf%phosphorus%ptr
    veg_litterfall_t%leaf%phosphorus  = veg_litterfall_t%leaf%phosphorus  - retranslocation_leaf_p

    !> 3.0 fine root turnover to litter
    !! 
    fturn = lctlib%tau_fine_root * dt/one_day/one_year
    veg_litterfall_t%fine_root = veg_pool%fine_root * fturn

    !> 4.0 coarse root turnover to litter 
    !! 
    fturn = lctlib%tau_coarse_root * dt/one_day/one_year 
    veg_litterfall_t%coarse_root = veg_pool%coarse_root * fturn

    !> 5.0 sapwood turnover 
    !! 

    !> 5.1 litter production from halms (grasses), or branches (trees)
    !! 
    SELECT CASE (lctlib%growthform)
      CASE (itree)
        fturn = lctlib%tau_branches * dt/one_day/one_year 
      CASE (igrass) ! for grases follows leaf phenology
        fturn = fturn_leaf
    ENDSELECT 
    veg_litterfall_t%sap_wood = veg_pool%sap_wood * fturn

    !> 5.2 heart wood formation, including nutrient recycling
    !! 
    IF(lctlib%growthform == itree)THEN
        fturn = lctlib%tau_sap_wood * dt/one_day/one_year 
        veg_growth_t%heart_wood = veg_pool%sap_wood     * fturn
        veg_growth_t%sap_wood   = veg_growth_t%sap_wood - veg_growth_t%heart_wood 

        ! Nutrient recycling
        nitrogen_recycling_fract            = calc_nitrogen_recycling_fract()
        phosphorus_recycling_fract          = calc_phosphorus_recycling_fract()
        retranslocation_sap_wood_n          = veg_growth_t%heart_wood%nitrogen%ptr   * nitrogen_recycling_fract
        veg_growth_t%heart_wood%nitrogen    = veg_growth_t%heart_wood%nitrogen       * (1._wp - nitrogen_recycling_fract)
        retranslocation_sap_wood_n15        = veg_growth_t%heart_wood%nitrogen15%ptr * nitrogen_recycling_fract
        veg_growth_t%heart_wood%nitrogen15  = veg_growth_t%heart_wood%nitrogen15     * (1._wp - nitrogen_recycling_fract)
        retranslocation_sap_wood_p          = veg_growth_t%heart_wood%phosphorus%ptr * phosphorus_recycling_fract
        veg_growth_t%heart_wood%phosphorus  = veg_growth_t%heart_wood%phosphorus     * (1._wp - phosphorus_recycling_fract)
    ENDIF 

    !> 6.0 fruit turnover to seed_bed
    !! 
    fturn                 = lctlib%tau_fruit   * dt/one_day/one_year 
    veg_growth_t%seed_bed = veg_pool%fruit     * fturn
    veg_growth_t%fruit    = veg_growth_t%fruit - veg_growth_t%seed_bed

    !> 7.0 seed-bed turnover to litter 
    !! 
    fturn                     = lctlib%tau_seed_litter * dt/one_day/one_year 
    veg_litterfall_t%seed_bed = veg_pool%seed_bed      * fturn

    !> 8.0 update litter, growth and recycling of the plant class
    !!
    veg_litterfall             = veg_litterfall + veg_litterfall_t
    veg_growth                 = veg_growth     + veg_growth_t

    recycling_leaf_n       = recycling_leaf_n       + retranslocation_leaf_n
    recycling_leaf_n15     = recycling_leaf_n15     + retranslocation_leaf_n15
    recycling_leaf_p       = recycling_leaf_p       + retranslocation_leaf_p
    recycling_sap_wood_n   = recycling_sap_wood_n   + retranslocation_sap_wood_n
    recycling_sap_wood_n15 = recycling_sap_wood_n15 + retranslocation_sap_wood_n15
    recycling_sap_wood_p   = recycling_sap_wood_p   + retranslocation_sap_wood_p

  END SUBROUTINE calc_turnover


  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional N recycling from leaves 
  !!
  !! should be made dynamically dependended on N stress later 
  !!
  !! @todo veg_turnover:calc_nitrogen_recycling_fract implement functionality other than "set to default"
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_nitrogen_recycling_fract() RESULT (recycling)

    USE mo_veg_constants,     ONLY: resorption_default

    IMPLICIT NONE
    ! 0.1 InOut

    ! 0.2 Local
    REAL(wp) :: recycling


    recycling = resorption_default

  END FUNCTION calc_nitrogen_recycling_fract


  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional P recycling from leaves
  !!
  !! should be made dynamically dependended on P stress later
  !!
  !! @todo veg_turnover:calc_phosphorus_recycling_fract implement functionality other than "set to default"
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_phosphorus_recycling_fract() RESULT (recycling)

    USE mo_veg_constants,     ONLY: resorption_default

    IMPLICIT NONE
    ! 0.1 InOut

    ! 0.2 Local
    REAL(wp) :: recycling 


    recycling = resorption_default

  END FUNCTION calc_phosphorus_recycling_fract

END MODULE mo_veg_turnover
