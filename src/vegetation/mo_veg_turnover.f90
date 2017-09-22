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

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_turnover'

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
  SUBROUTINE update_veg_turnover(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg        !  USE mo_veg_memory_class,      ONLY: t_veg_element_a, t_veg_pool_a
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)    :: options          !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_turnover'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                        ! TYPE(t_pheno_memory),    POINTER :: pheno__mem
    ! Declare pointers to variables in memory
    ! PHENO_
    dsl4jsb_Logical2D_onChunk   :: growing_season
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: target_lai                      ! VEG_ 
    dsl4jsb_Real2D_onChunk      :: target_cn_fine_root
    dsl4jsb_Real2D_onChunk      :: target_np_fine_root
    dsl4jsb_Real2D_onChunk      :: target_cn_leaf
    dsl4jsb_Real2D_onChunk      :: target_np_leaf
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n15
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    ! Declare pointers to structures in memory
#define Def_veg_pool
#define Def_veg_litterfall
#define Def_veg_growth
#include "slm_pools_a_incl.f90"
#undef Def_veg_pool
#undef Def_veg_litterfall
#undef Def_veg_growth
!    TYPE(t_veg_memory),       POINTER     :: veg      ! local pointer to veg__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)                        ! veg__mem => veg_memory
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)           ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_fine_root) 
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_fine_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)
    ! Get process memory structures
#define Get_veg_pool
#define Get_veg_litterfall
#define Get_veg_growth
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#undef Get_veg_litterfall
#undef Get_veg_growth

!! *slmtemp
!    veg     => veg__mem 



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Calculate nutrient flux from enzyme, RNA and other turnover in leaves and roots
    !!
    CALL calc_nutrient_recycling( &
                                  veg_pool                              , &
                                  target_cn_fine_root(:)                , &
                                  target_np_fine_root(:)                , &
                                  target_cn_leaf(:)                     , &
                                  target_np_leaf(:)                     , &
                                  growing_season(:)                     , &
                                  recycling_leaf_n(:)                   , & ! inout
                                  recycling_fine_root_n(:)              , &
                                  recycling_leaf_p(:)                   , &
                                  recycling_fine_root_p(:)              , &
                                  recycling_leaf_n15(:)                 , &
                                  recycling_fine_root_n15(:)            )

    !> 2.0 Calculate recycling and retranslocation
    !!  
    CALL calc_turnover(lctlib                                       , & ! in 
                              veg_pool                              , & ! inout
                              veg_growth                            , & ! inout
                              veg_litterfall                        , & ! inout
                              growing_season(1)                     , & ! in
                              lai(1)                                , & ! in
                              target_lai(1)                         , & ! in
                              recycling_leaf_n(1)                   , & ! inout
                              recycling_leaf_n15(1)                 , &
                              recycling_leaf_p(1)                   , &
                              recycling_heart_wood_n(1)             , &
                              recycling_heart_wood_n15(1)           , &
                              recycling_heart_wood_p(1)             )


#define Dealloc_veg_pool
#define Dealloc_veg_litterfall
#define Dealloc_veg_growth
#include "slm_pools_a_incl.f90"
#undef Dealloc_veg_pool
#undef Dealloc_veg_litterfall
#undef Dealloc_veg_growth

  END SUBROUTINE update_veg_turnover


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_turnover
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Calculate nutrient flux from enzyme turnover and other biochemical turnover in leaves and roots
  !! at the time scale of tau_nutrient_recycling
  !!
  !! Calculations only during growing season \n
  !! Input: plant status and leaf and root C:N:P target ratios \n
  !! Output: net flux between labile and leaves/roots (+ve: flux to labile, -ve: flux from labile)
  !----------------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE calc_nutrient_recycling( &
                                                    veg_pool , &
                                                    target_cn_fine_root               , &
                                                    target_np_fine_root               , &
                                                    target_cn_leaf                    , &
                                                    target_np_leaf                    , &
                                                    growing_season                    , &
                                                    recycling_leaf_n                  , &
                                                    recycling_fine_root_n             , &
                                                    recycling_leaf_p                  , &
                                                    recycling_fine_root_p             , &
                                                    recycling_leaf_n15                , &
                                                    recycling_fine_root_n15           )

    USE mo_base_constants,        ONLY: dtime, one_day
    USE mo_veg_constants,         ONLY: tau_nutrient_recycling, tau_labile
    dsl4slm_Use_memory_a_veg


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_veg_pool_a), INTENT(in)    :: veg_pool
    REAL(wp),           INTENT(in)    :: target_cn_fine_root               , &  !< target C:N of fine_roots
                                         target_np_fine_root               , &  !< target N:P of fine_roots 
                                         target_cn_leaf                    , &  !< ..
                                         target_np_leaf                         !< ..
    LOGICAL,            INTENT(in)    :: growing_season                         !< whether the plant is in the growing season
    REAL(wp),           INTENT(inout) :: recycling_leaf_n                  , & 
                                         recycling_fine_root_n             , & 
                                         recycling_leaf_p                  , & 
                                         recycling_fine_root_p             , & 
                                         recycling_leaf_n15                , & 
                                         recycling_fine_root_n15
    ! 0.2 Local
    REAL(wp)                          :: labile_dn15              ! labile pools 15N to N molar mixing ratio
    REAL(wp)                          :: scal                     ! helper variable
    REAL(wp)                          :: recycling_leaf, &        ! more helpers1
                                         recycling_fine_root      ! more helpers2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nutrient_recycling'


    IF(growing_season)THEN 
       !> 1.0 replace overturning N from leaves and roots by new N from the labile pool using the target C:N
       !!
       recycling_leaf      = ( veg_pool%leaf%nitrogen - & 
                               veg_pool%leaf%carbon / target_cn_leaf ) / tau_nutrient_recycling/one_day * dtime
       recycling_fine_root = ( veg_pool%fine_root%nitrogen - &
                               veg_pool%fine_root%carbon / target_cn_fine_root ) / tau_nutrient_recycling/one_day * dtime
       !> 1.1 if N is addded to leaves or fine root (negative flux), make sure there is enough labile N to support this
       !! 
       IF( -1.0_wp * (recycling_leaf + recycling_fine_root) > eps8) THEN
          scal = -1.0_wp * ( 1.0_wp/tau_labile/one_day * veg_pool%labile%nitrogen * dtime ) / & 
                           ( recycling_leaf + recycling_fine_root ) 
          IF(scal < 1.0_wp) THEN
             recycling_leaf       = recycling_leaf * scal
             recycling_fine_root  = recycling_fine_root * scal
          ENDIF
       ENDIF
       recycling_leaf_n      = recycling_leaf_n      + recycling_leaf
       recycling_fine_root_n = recycling_fine_root_n + recycling_fine_root

       !> 1.2 calculate signature of N transfer from labile to leaf or inverse
       !!  
       IF(recycling_leaf >=0.0_wp) THEN
          IF(veg_pool%leaf%nitrogen > eps8) THEN 
             recycling_leaf_n15 = recycling_leaf_n15 + & 
               recycling_leaf/veg_pool%leaf%nitrogen * veg_pool%leaf%nitrogen15
          ENDIF
          IF(veg_pool%labile%nitrogen > eps8) THEN 
             recycling_leaf_n15 = recycling_leaf_n15 + & 
               recycling_leaf/veg_pool%labile%nitrogen * veg_pool%labile%nitrogen15
          ENDIF
       ENDIF
       IF(recycling_fine_root >=0.0_wp) THEN
          IF(veg_pool%fine_root%nitrogen > eps8) THEN 
             recycling_fine_root_n15 = recycling_fine_root_n15 + & 
                recycling_fine_root/veg_pool%fine_root%nitrogen * veg_pool%fine_root%nitrogen15
          ENDIF
          IF(veg_pool%labile%nitrogen > eps8) THEN 
             recycling_fine_root_n15 = recycling_fine_root_n15 + & 
                recycling_fine_root/veg_pool%labile%nitrogen * veg_pool%labile%nitrogen15
          ENDIF
       ENDIF 

       !> 2.0 replace overturning P from leaves and roots by new N from the labile pool using the target N:P
       !!
       recycling_leaf      = ( veg_pool%leaf%phosphorus - & 
                               veg_pool%leaf%nitrogen / target_np_leaf ) / tau_nutrient_recycling/one_day * dtime 
       recycling_fine_root = ( veg_pool%fine_root%phosphorus - &
                               veg_pool%fine_root%nitrogen / target_np_fine_root ) / tau_nutrient_recycling/one_day * dtime
       !> 2.1 if P is addded to leaves or fine root, make sure there is enough labile P to support this
       !! 
       IF( -1.0_wp * ( recycling_leaf + recycling_fine_root ) > eps8) THEN
          scal = -1.0_wp * ( 1.0_wp/tau_labile/one_day * veg_pool%labile%phosphorus * dtime ) / & 
                           ( recycling_leaf + recycling_fine_root ) 
          IF(scal < 1.0_wp) THEN
             recycling_leaf       = recycling_leaf * scal
             recycling_fine_root  = recycling_fine_root * scal
          ENDIF
       ENDIF
       recycling_leaf_p        = recycling_leaf_p        + recycling_leaf
       recycling_fine_root_p   = recycling_fine_root_p   + recycling_fine_root
   
    ENDIF

  END SUBROUTINE calc_nutrient_recycling


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_turnover
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover: recycling and retranslocation 
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
!  SUBROUTINE calc_turnover(lctlib                     , &
  ELEMENTAL SUBROUTINE calc_turnover(lctlib                     , &
                                      veg_pool                  , &
                                      veg_growth                , &
                                      veg_litterfall            , &
                                      growing_season            , &
                                      lai                       , &
                                      target_lai                , &
                                      recycling_leaf_n          , &
                                      recycling_leaf_n15        , &
                                      recycling_leaf_p          , &
                                      recycling_heart_wood_n    , &
                                      recycling_heart_wood_n15  , &
                                      recycling_heart_wood_p      )

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_base_constants,        ONLY: dtime, one_day, one_year
    USE mo_shared_constants
    USE mo_veg_constants
    USE mo_pheno_constants,       ONLY: ievergreen, iraingreen, isummergreen, iperennial
    dsl4slm_Use_memory_a_veg

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_lctlib_element),       INTENT(in)    :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_veg_pool_a),           INTENT(inout) :: veg_pool         , & !< current state of vegetation pools 
                                                   veg_growth       , & !< current growth rate of vegetation pools
                                                   veg_litterfall       !< current litter fall from vegetation pools
    LOGICAL,                      INTENT(in)    :: growing_season       !< whether the plant is in the growing season
    REAL(wp),                     INTENT(in)    :: lai                  !< the plant's LAI
    REAL(wp),                     INTENT(in)    :: target_lai           !< the plant's annual maximum LAI
    REAL(wp),                     INTENT(inout) :: recycling_leaf_n , & !< foliar N returning to labile pool 
                                                   recycling_leaf_n15, & !< foliar N15 returning to labile pool 
                                                   recycling_leaf_p , & !< foliar P returning to labile pool 
                                                   recycling_heart_wood_n    , & !< woody N returning to labile pool 
                                                   recycling_heart_wood_n15  , & !< woody N15 returning to labile pool
                                                   recycling_heart_wood_p        !< woody P returning to labile pool
    ! 0.2 Local
    REAL(wp)                          :: retranslocation              , & !< retranslocation before shedding 
                                         fturn_leaf                   , & !< fraction of leaves turning over 
                                         fturn                            !< other tissue turnover fraction
    REAL(wp)                          :: nitrogen_recycling_fract     , & !< recycling fraction of N
                                         phosphorus_recycling_fract       !< recycling fraction of P
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_turnover'


    !> 1.0 foliage turnover, including recycling of nutrients 
    !!
    !! flat rate for evergreens, dual rates for deciduous (senescence after the growing season)
    SELECT CASE (lctlib%phenology_type)
      CASE (ievergreen)
         fturn_leaf = 1.0_wp/lctlib%tau_leaf/one_day/one_year * dtime
      CASE (iraingreen, isummergreen)
         IF(growing_season)THEN
            fturn_leaf = 0.0_wp              ! herbivory
         ELSE
             IF(lai>eps4)THEN
                fturn_leaf = MAX(eps4,max_leaf_shedding_rate*dtime/one_day*target_lai/lai)
             ELSE
                fturn_leaf = 1.0_wp
             ENDIF
         ENDIF
      CASE (iperennial)
         IF(growing_season)THEN
            fturn_leaf = 1.0_wp/lctlib%tau_leaf/one_day/one_year * dtime
         ELSE
             IF(lai>eps4)THEN
                fturn_leaf = MAX(1.0_wp/lctlib%tau_leaf/one_day/one_year * dtime, & 
                                 max_leaf_shedding_rate*dtime/one_day*target_lai/lai)
             ELSE
                fturn_leaf = 1.0_wp
             ENDIF
       ENDIF
    END SELECT
    veg_litterfall%leaf = veg_litterfall%leaf + veg_pool%leaf * fturn_leaf

    !> 1.1 Nutrient recycling from foliage turnover
    !! 
    nitrogen_recycling_fract    = calc_nitrogen_recycling_fract("leaf")
    phosphorus_recycling_fract  = calc_phosphorus_recycling_fract("leaf")

    retranslocation = nitrogen_recycling_fract * veg_pool%leaf%nitrogen * fturn_leaf 
    recycling_leaf_n = recycling_leaf_n + retranslocation
    veg_litterfall%leaf%nitrogen = veg_litterfall%leaf%nitrogen - retranslocation

    retranslocation = nitrogen_recycling_fract * veg_pool%leaf%nitrogen15 * fturn_leaf
    recycling_leaf_n15 = recycling_leaf_n15 + retranslocation
    veg_litterfall%leaf%nitrogen15 = veg_litterfall%leaf%nitrogen15 - retranslocation

    retranslocation = phosphorus_recycling_fract * veg_pool%leaf%phosphorus * fturn_leaf 
    recycling_leaf_p = recycling_leaf_p + retranslocation
    veg_litterfall%leaf%phosphorus = veg_litterfall%leaf%phosphorus - retranslocation

    !> 2.0 fine root turnover to litter
    !! 
    fturn = 1.0_wp/lctlib%tau_fine_root/one_day/one_year * dtime
    veg_litterfall%fine_root = veg_litterfall%fine_root + veg_pool%fine_root * fturn

    !> 3.0 coarse root turnover to litter 
    !! 
    fturn = 1.0_wp/lctlib%tau_coarse_root/one_day/one_year * dtime
    veg_litterfall%coarse_root = veg_litterfall%coarse_root + veg_pool%coarse_root * fturn

    !> 4.0 sapwood turnover 
    !! 

    !> 4.1 litter production from halms (grasses), or branches (trees)
    !! 
    SELECT CASE (lctlib%growthform)
      CASE (itree)
        fturn = 1.0_wp/lctlib%tau_branch/one_day/one_year * dtime
      CASE (igrass) ! for grases follows leaf phenology
        fturn = fturn_leaf
    ENDSELECT 
    veg_litterfall%sap_wood = veg_litterfall%sap_wood + lctlib%frac_sapwood_branch * veg_pool%sap_wood * fturn

    !> 4.2 heart wood formation, including nutrient recycling
    !! 
    IF(lctlib%growthform == itree)THEN
        fturn = 1.0_wp/lctlib%tau_sap_wood/one_day/one_year * dtime
        veg_growth%heart_wood = veg_growth%heart_wood + veg_pool%sap_wood * fturn
        veg_growth%sap_wood   = veg_growth%sap_wood   - veg_pool%sap_wood * fturn 

        ! Nutrient recycling
        nitrogen_recycling_fract            = calc_nitrogen_recycling_fract("wood")
        phosphorus_recycling_fract          = calc_phosphorus_recycling_fract("wood")
        recycling_heart_wood_n   = recycling_heart_wood_n   + & 
                                   veg_pool%sap_wood%nitrogen   * fturn * nitrogen_recycling_fract
        recycling_heart_wood_n15 = recycling_heart_wood_n15 + &
                                   veg_pool%sap_wood%nitrogen15 * fturn * nitrogen_recycling_fract
        recycling_heart_wood_p   = recycling_heart_wood_p   + & 
                                   veg_pool%sap_wood%phosphorus * fturn * phosphorus_recycling_fract
    ENDIF 

    !> 5.0 fruit turnover to seed_bed
    !! 
    fturn               = 1.0_wp/lctlib%tau_fruit/one_day/one_year * dtime
    veg_growth%seed_bed = veg_growth%seed_bed + veg_pool%fruit * fturn
    veg_growth%fruit    = veg_growth%fruit    - veg_pool%fruit * fturn 

    !> 6.0 seed-bed turnover to litter 
    !! 
    fturn                   = 1.0_wp/lctlib%tau_seed_litter/one_day/one_year * dtime
    veg_litterfall%seed_bed = veg_litterfall%seed_bed + veg_pool%seed_bed    * fturn

  END SUBROUTINE calc_turnover


  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional N recycling from leaves 
  !!
  !! should be made dynamically dependended on N stress later 
  !!
  !! @todo veg_turnover:calc_nitrogen_recycling_fract implement functionality other than "set to default"
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_nitrogen_recycling_fract(pool_name) RESULT (recycling_fract)

    USE mo_veg_constants,     ONLY: resorp_fract_leaf,resorp_fract_wood

    IMPLICIT NONE
    ! 0.1 InOut
    CHARACTER(4), INTENT(in) :: pool_name
    ! 0.2 Local
    REAL(wp) :: recycling_fract 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nitrogen_recycling_fract'

    SELECT CASE(TRIM(pool_name))
      CASE("leaf")
        recycling_fract = resorp_fract_leaf 
      CASE("wood")
        recycling_fract = resorp_fract_wood
    ENDSELECT


  END FUNCTION calc_nitrogen_recycling_fract


  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional P recycling from leaves
  !!
  !! should be made dynamically dependended on P stress later
  !!
  !! @todo veg_turnover:calc_phosphorus_recycling_fract implement functionality other than "set to default"
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_phosphorus_recycling_fract(pool_name) RESULT (recycling_fract)

    USE mo_veg_constants,     ONLY: resorp_fract_leaf,resorp_fract_wood

    IMPLICIT NONE
    ! 0.1 InOut
    CHARACTER(4), INTENT(in) :: pool_name
    ! 0.2 Local
    REAL(wp) :: recycling_fract 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_phosphorus_recycling_fract'

    SELECT CASE(TRIM(pool_name))
      CASE("leaf")
        recycling_fract = resorp_fract_leaf 
      CASE("wood")
        recycling_fract = resorp_fract_wood
    ENDSELECT

  END FUNCTION calc_phosphorus_recycling_fract

END MODULE mo_veg_turnover
