!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
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

#ifdef __QUINCY_NO_JSB__
#include "dsl_jsbach.inc"
#endif

  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message

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

    USE mo_kind,                  ONLY: wp
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    ! Use of pools
#ifdef __QUINCY_NO_JSB__
    dsl4slm_Use_memory_a_veg        !  USE mo_veg_memory_class,      ONLY: t_veg_element_a, t_veg_pool_a
    USE mo_veg_memory_class,      ONLY: veg_pool, veg_growth, veg_litterfall
#else
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),        POINTER       :: model
    TYPE(t_lctlib_element),   POINTER       :: lctlib           !< land-cover-type library - parameter across pft's
    INTEGER                                 :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_turnover'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration, process memory and pools
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                        ! TYPE(t_pheno_memory),    POINTER :: pheno__mem
#ifdef __QUINCY_NO_JSB__
#else
    dsl4jsb_Def_pool(t_veg_tissues_pool)  :: veg_pool
    dsl4jsb_Def_pool(t_veg_tissues_pool)  :: veg_growth
    dsl4jsb_Def_pool(t_veg_tissues_pool)  :: veg_litterfall
#endif
    ! Declare pointers to variables in memory
    ! PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season
    dsl4jsb_Real2D_onChunk      :: lai
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: target_lai
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
    dsl4jsb_Real2D_onChunk      :: net_growth 
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! Get process config
!    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)                        ! veg__mem => veg_memory
    dsl4jsb_Get_memory(PHENO_)
#ifdef __QUINCY_NO_JSB__
#else
    ! Get pools
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_pool,          veg_pool)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_growth,        veg_growth)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_litterfall,    veg_litterfall)
    ! associate element-variables with chunk
    CALL veg_pool%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_growth%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_litterfall%Associate_chunk_var_ptr(ics, ice, iblk)
#endif
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)             ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                        ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, target_lai)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_fine_root)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_fine_root)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_leaf)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_leaf)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)           ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)      ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)       ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)     ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)       ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth)                   ! out
#ifdef __QUINCY_NO_JSB__
    ! Get pools
#define Get_veg_pool
#define Get_veg_litterfall
#define Get_veg_growth
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#undef Get_veg_litterfall
#undef Get_veg_growth
#endif


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Calculate nutrient flux from enzyme, RNA and other turnover in leaves and roots
    !!
    CALL calc_nutrient_recycling( nc                                    , & ! in
                                  target_cn_fine_root(:)                , &
                                  target_np_fine_root(:)                , &
                                  target_cn_leaf(:)                     , &
                                  target_np_leaf(:)                     , &
                                  growing_season(:)                     , &
                                  veg_pool                              , & ! in
                                  recycling_leaf_n(:)                   , & ! inout
                                  recycling_fine_root_n(:)              , &
                                  recycling_leaf_p(:)                   , &
                                  recycling_fine_root_p(:)              , &
                                  recycling_leaf_n15(:)                 , &
                                  recycling_fine_root_n15(:)            )   ! inout

    !> 2.0 Calculate recycling and retranslocation
    !!  
    CALL calc_turnover( nc                                    , & ! in 
                        lctlib%phenology_type                 , &
                        lctlib%growthform                     , &
                        lctlib%frac_sapwood_branch            , &
                        lctlib%tau_leaf                       , &
                        lctlib%tau_fine_root                  , &
                        lctlib%tau_coarse_root                , &
                        lctlib%tau_branch                     , &
                        lctlib%tau_sap_wood                   , &
                        lctlib%tau_fruit                      , &
                        lctlib%tau_seed_litter                , &
                        growing_season(:)                     , & ! in
                        lai(:)                                , & ! in
                        target_lai(:)                         , & ! in
                        recycling_leaf_n(:)                   , & ! inout
                        recycling_leaf_n15(:)                 , &
                        recycling_leaf_p(:)                   , &
                        recycling_heart_wood_n(:)             , &
                        recycling_heart_wood_n15(:)           , &
                        recycling_heart_wood_p(:)             , &
                        veg_pool                              , &  ! in
                        veg_growth                            , &  ! inout
                        veg_litterfall                        , &  ! inout
                        net_growth                             )   ! out

#ifdef __QUINCY_NO_JSB__
! Pass back values to process memory
#define Passback_veg_pool
#define Passback_veg_litterfall
#define Passback_veg_growth
#include "slm_pools_a_incl.f90"
#undef Passback_veg_pool
#undef Passback_veg_litterfall
#undef Passback_veg_growth
#endif

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
  SUBROUTINE calc_nutrient_recycling( nc                                , &
                                      target_cn_fine_root               , &
                                      target_np_fine_root               , &
                                      target_cn_leaf                    , &
                                      target_np_leaf                    , &
                                      growing_season                    , &
                                      veg_pool                          , &
                                      recycling_leaf_n                  , &
                                      recycling_fine_root_n             , &
                                      recycling_leaf_p                  , &
                                      recycling_fine_root_p             , &
                                      recycling_leaf_n15                , &
                                      recycling_fine_root_n15           )

    USE mo_kind,                  ONLY: wp
    USE mo_jsb_impl_constants,    ONLY: test_false_true
    USE mo_jsb_math_constants,    ONLY: dtime, one_day, eps8
    USE mo_veg_constants,         ONLY: tau_nutrient_recycling, tau_labile
#ifdef __QUINCY_NO_JSB__
    dsl4slm_Use_memory_a_veg
#else
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                          INTENT(in)    :: nc                           !< dimensions
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: target_cn_fine_root               , &  !< target C:N of fine_roots
                                                       target_np_fine_root               , &  !< target N:P of fine_roots 
                                                       target_cn_leaf                    , &  !< ..
                                                       target_np_leaf                         !< ..
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: growing_season                         !< whether the plant is in the growing season
#ifdef __QUINCY_NO_JSB__
    TYPE(t_veg_pool_a),               INTENT(in)    :: veg_pool                               !< vegetation pool
#else
    TYPE(t_veg_tissues_pool),         INTENT(in)    :: veg_pool                               !< vegetation pool
#endif
    REAL(wp), DIMENSION(nc),          INTENT(inout) :: recycling_leaf_n                  , & 
                                                       recycling_fine_root_n             , & 
                                                       recycling_leaf_p                  , & 
                                                       recycling_fine_root_p             , & 
                                                       recycling_leaf_n15                , & 
                                                       recycling_fine_root_n15
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                         :: hlp1                     ! helper variable
    REAL(wp), DIMENSION(nc)                         :: recycling_leaf, &        ! more helpers1
                                                       recycling_fine_root      ! more helpers2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nutrient_recycling'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 replace overturning N from leaves and roots by new N from the labile pool using the target C:N
    !!
    WHERE (growing_season(:) > test_false_true)
       recycling_leaf(:)      = (veg_pool%leaf%nitrogen(:) - & 
                                 veg_pool%leaf%carbon(:) / target_cn_leaf(:)) / tau_nutrient_recycling / &
                                 one_day * dtime
       recycling_fine_root(:) = (veg_pool%fine_root%nitrogen(:) - &
                                 veg_pool%fine_root%carbon(:) / target_cn_fine_root(:)) / tau_nutrient_recycling / &
                                 one_day * dtime
       !> 1.1 if N is addded to leaves or fine root (negative flux), make sure there is enough labile N to support this
       !! 
       WHERE((-1.0_wp * (recycling_leaf(:) + recycling_fine_root(:))) > eps8)
          hlp1(:) = -1.0_wp * (1.0_wp / tau_labile / one_day * veg_pool%labile%nitrogen(:) * dtime) / & 
                              (recycling_leaf(:) + recycling_fine_root(:)) 
          WHERE(hlp1(:) < 1.0_wp)
             recycling_leaf(:)       = recycling_leaf(:)      * hlp1(:)
             recycling_fine_root(:)  = recycling_fine_root(:) * hlp1(:)
          END WHERE
       END WHERE
       recycling_leaf_n(:)      = recycling_leaf_n(:)      + recycling_leaf(:)
       recycling_fine_root_n(:) = recycling_fine_root_n(:) + recycling_fine_root(:)

       !> 1.2 calculate signature of N transfer from labile to leaf or inverse
       !!  
       WHERE(recycling_leaf(:) >= 0.0_wp)  ! for all positive values of recycling_fine_root(:)
          WHERE(veg_pool%leaf%nitrogen(:) > eps8)
            recycling_leaf_n15(:) = recycling_leaf_n15(:) + & 
                                    recycling_leaf(:) / veg_pool%leaf%nitrogen(:) * veg_pool%leaf%nitrogen15(:)
          END WHERE
       ELSEWHERE  ! for all negative values of recycling_fine_root(:)
          WHERE(veg_pool%labile%nitrogen(:) > eps8)
            recycling_leaf_n15(:) = recycling_leaf_n15(:) + & 
                                    recycling_leaf(:) / veg_pool%labile%nitrogen(:) * veg_pool%labile%nitrogen15(:)
          END WHERE
       END WHERE
       WHERE(recycling_fine_root(:) >= 0.0_wp)  ! for all positive values of recycling_fine_root(:)
          WHERE(veg_pool%fine_root%nitrogen(:) > eps8)
            recycling_fine_root_n15(:) = recycling_fine_root_n15(:) + & 
                                         recycling_fine_root(:) / veg_pool%fine_root%nitrogen(:) * veg_pool%fine_root%nitrogen15(:)
          END WHERE
       ELSEWHERE  ! for all negative values of recycling_fine_root(:)
          WHERE(veg_pool%labile%nitrogen(:) > eps8)
            recycling_fine_root_n15(:) = recycling_fine_root_n15(:) + & 
                                         recycling_fine_root(:) / veg_pool%labile%nitrogen(:) * veg_pool%labile%nitrogen15(:)
          END WHERE
       END WHERE


       !> 2.0 replace overturning P from leaves and roots by new N from the labile pool using the target N:P
       !!
       recycling_leaf(:)      = (veg_pool%leaf%phosphorus(:) - veg_pool%leaf%nitrogen(:) / target_np_leaf(:)) / &
                                 tau_nutrient_recycling / one_day * dtime 
       recycling_fine_root(:) = (veg_pool%fine_root%phosphorus(:) - veg_pool%fine_root%nitrogen(:) / target_np_fine_root(:)) / &
                                 tau_nutrient_recycling / one_day * dtime
       !> 2.1 if P is addded to leaves or fine root, make sure there is enough labile P to support this
       !! 
       WHERE((-1.0_wp * (recycling_leaf(:) + recycling_fine_root(:))) > eps8)
          hlp1(:) = -1.0_wp * (1.0_wp/tau_labile/one_day * veg_pool%labile%phosphorus(:) * dtime) / & 
                              (recycling_leaf(:) + recycling_fine_root(:)) 
          WHERE(hlp1(:) < 1.0_wp)
             recycling_leaf(:)       = recycling_leaf(:)      * hlp1(:)
             recycling_fine_root(:)  = recycling_fine_root(:) * hlp1(:)
          END WHERE
       END WHERE
       recycling_leaf_p(:)        = recycling_leaf_p(:)        + recycling_leaf(:)
       recycling_fine_root_p(:)   = recycling_fine_root_p(:)   + recycling_fine_root(:)
   
    END WHERE ! growing_season(:) > test_false_true

  END SUBROUTINE calc_nutrient_recycling


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_turnover
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover: recycling and retranslocation 
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_turnover( nc                        , &
                            lctlib_phenology_type     , &
                            lctlib_growthform         , &
                            lctlib_frac_sapwood_branch , &
                            lctlib_tau_leaf           , &
                            lctlib_tau_fine_root      , &
                            lctlib_tau_coarse_root    , &
                            lctlib_tau_branch         , &
                            lctlib_tau_sap_wood       , &
                            lctlib_tau_fruit          , &
                            lctlib_tau_seed_litter    , &
                            growing_season            , &
                            lai                       , &
                            target_lai                , &
                            recycling_leaf_n          , &
                            recycling_leaf_n15        , &
                            recycling_leaf_p          , &
                            recycling_heart_wood_n    , &
                            recycling_heart_wood_n15  , &
                            recycling_heart_wood_p    , &  
                            veg_pool                  , &
                            veg_growth                , &
                            veg_litterfall            , &
                            net_growth                )

    USE mo_kind,                  ONLY: wp
    USE mo_jsb_impl_constants,    ONLY: test_false_true
    USE mo_jsb_math_constants,    ONLY: dtime, one_day, one_year, eps4
    USE mo_veg_constants,         ONLY: max_leaf_shedding_rate, igrass, itree
    USE mo_pheno_constants,       ONLY: ievergreen, iraingreen, isummergreen, iperennial
#ifdef __QUINCY_NO_JSB__
    dsl4slm_Use_memory_a_veg
#else
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                         INTENT(in)    :: nc                           !< dimensions
    INTEGER,                         INTENT(in)    :: lctlib_phenology_type        !< lctlib parameter
    INTEGER,                         INTENT(in)    :: lctlib_growthform            !< lctlib parameter
    REAL(wp),                        INTENT(in)    :: lctlib_frac_sapwood_branch   !< lctlib parameter
    REAL(wp),                        INTENT(in)    :: lctlib_tau_leaf              !< lctlib parameter
    REAL(wp),                        INTENT(in)    :: lctlib_tau_fine_root         !< lctlib parameter
    REAL(wp),                        INTENT(in)    :: lctlib_tau_coarse_root       !< lctlib parameter
    REAL(wp),                        INTENT(in)    :: lctlib_tau_branch            !< lctlib parameter
    REAL(wp),                        INTENT(in)    :: lctlib_tau_sap_wood          !< lctlib parameter
    REAL(wp),                        INTENT(in)    :: lctlib_tau_fruit             !< lctlib parameter
    REAL(wp),                        INTENT(in)    :: lctlib_tau_seed_litter       !< lctlib parameter
    REAL(wp), DIMENSION(nc),         INTENT(in)    :: growing_season               !< whether the plant is in the growing season
    REAL(wp), DIMENSION(nc),         INTENT(in)    :: lai                      , & !< the plant's LAI
                                                      target_lai                   !< the plant's annual maximum LAI
    REAL(wp), DIMENSION(nc),         INTENT(inout) :: recycling_leaf_n         , & !< foliar N returning to labile pool 
                                                      recycling_leaf_n15       , & !< foliar N15 returning to labile pool 
                                                      recycling_leaf_p         , & !< foliar P returning to labile pool 
                                                      recycling_heart_wood_n   , & !< woody N returning to labile pool 
                                                      recycling_heart_wood_n15 , & !< woody N15 returning to labile pool
                                                      recycling_heart_wood_p       !< woody P returning to labile pool
#ifdef __QUINCY_NO_JSB__
    TYPE(t_veg_pool_a),              INTENT(in)    :: veg_pool                     !< current state of vegetation pools 
    TYPE(t_veg_pool_a),              INTENT(inout) :: veg_growth                   !< current growth rate of vegetation pools
    TYPE(t_veg_pool_a),              INTENT(inout) :: veg_litterfall               !< current litter fall from vegetation pools
#else
    TYPE(t_veg_tissues_pool),        INTENT(in)    :: veg_pool                     !< current state of vegetation pools 
    TYPE(t_veg_tissues_pool),        INTENT(inout) :: veg_growth                   !< current growth rate of vegetation pools
    TYPE(t_veg_tissues_pool),        INTENT(inout) :: veg_litterfall               !< current litter fall from vegetation pools
#endif
    REAL(wp), DIMENSION(nc),         INTENT(out)   :: net_growth                   !< current net growth (µmol C / m2 / s) 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)           :: fturn_leaf                   , & !< fraction of leaves turning over 
                                         fturn                        , & !< other tissue turnover fraction
                                         retranslocation              , & !< retranslocation before shedding 
                                         nitrogen_recycling_fract     , & !< recycling fraction of N
                                         phosphorus_recycling_fract       !< recycling fraction of P
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_turnover'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !>0.9 init output variable
    !>
    net_growth(:) = 0.0_wp

    !> 1.0 foliage turnover, including recycling of nutrients 
    !!
    !! flat rate for evergreens, dual rates for deciduous (senescence after the growing season)
    SELECT CASE (lctlib_phenology_type)
      CASE (ievergreen)
         fturn_leaf(:)      = 1.0_wp / lctlib_tau_leaf / one_day / one_year * dtime
      CASE (iraingreen, isummergreen)
         WHERE(growing_season(:) > test_false_true)
            fturn_leaf(:)   = 0.0_wp              ! herbivory
         ELSEWHERE
           WHERE(lai(:) > eps4)
              fturn_leaf(:) = MAX(eps4, max_leaf_shedding_rate * dtime / one_day * target_lai(:) / lai(:))
           ELSEWHERE
              fturn_leaf(:) = 1.0_wp
           END WHERE
         END WHERE
      CASE (iperennial)
         WHERE(growing_season(:) > test_false_true)
           fturn_leaf(:)    = 1.0_wp / lctlib_tau_leaf / one_day / one_year * dtime
         ELSEWHERE
           WHERE(lai(:) > eps4)
             fturn_leaf(:)  = MAX(1.0_wp / lctlib_tau_leaf / one_day / one_year * dtime, & 
                                 max_leaf_shedding_rate * dtime / one_day * target_lai(:) / lai(:))
           ELSEWHERE
             fturn_leaf(:)  = 1.0_wp
           END WHERE
         END WHERE
    END SELECT
    veg_litterfall%leaf = veg_litterfall%leaf + veg_pool%leaf * fturn_leaf(:)
    net_growth(:)       = net_growth(:) - veg_pool%leaf%carbon(:) * fturn_leaf(:)

    !> 1.1 Nutrient recycling from foliage turnover
    !! 
    nitrogen_recycling_fract(:)       = calc_nitrogen_recycling_fract("leaf")
    phosphorus_recycling_fract(:)     = calc_phosphorus_recycling_fract("leaf")

    retranslocation(:)                = nitrogen_recycling_fract(:) * veg_pool%leaf%nitrogen(:) * fturn_leaf(:) 
    recycling_leaf_n(:)               = recycling_leaf_n(:) + retranslocation(:)
    veg_litterfall%leaf%nitrogen(:)   = veg_litterfall%leaf%nitrogen(:) - retranslocation(:)

    retranslocation(:)                = nitrogen_recycling_fract(:) * veg_pool%leaf%nitrogen15(:) * fturn_leaf(:)
    recycling_leaf_n15(:)             = recycling_leaf_n15(:) + retranslocation(:)
    veg_litterfall%leaf%nitrogen15(:) = veg_litterfall%leaf%nitrogen15(:) - retranslocation(:)

    retranslocation(:)                = phosphorus_recycling_fract(:) * veg_pool%leaf%phosphorus(:) * fturn_leaf(:) 
    recycling_leaf_p(:)               = recycling_leaf_p(:) + retranslocation(:)
    veg_litterfall%leaf%phosphorus(:) = veg_litterfall%leaf%phosphorus(:) - retranslocation(:)


    !> 2.0 fine root turnover to litter
    !! 
    fturn(:)                  = 1.0_wp / lctlib_tau_fine_root / one_day / one_year * dtime
    veg_litterfall%fine_root  = veg_litterfall%fine_root + veg_pool%fine_root * fturn(:)
    net_growth(:)             = net_growth(:) - veg_pool%fine_root%carbon(:) * fturn(:)


    !> 3.0 coarse root turnover to litter 
    !! 
    fturn(:)                    = 1.0_wp /lctlib_tau_coarse_root /one_day /one_year *dtime
    veg_litterfall%coarse_root  = veg_litterfall%coarse_root + veg_pool%coarse_root * fturn(:)
    net_growth(:)               = net_growth(:) - veg_pool%coarse_root%carbon(:) * fturn(:)


    !> 4.0 sapwood turnover 
    !! 
    !> 4.1 litter production from halms (grasses), or branches (trees)
    !! 
    SELECT CASE (lctlib_growthform)
      CASE (itree)
        fturn(:) = 1.0_wp / lctlib_tau_branch / one_day / one_year * dtime
      CASE (igrass) ! for grases follows leaf phenology
        fturn(:) = fturn_leaf(:)
    ENDSELECT 
    veg_litterfall%sap_wood = veg_litterfall%sap_wood + lctlib_frac_sapwood_branch * veg_pool%sap_wood * fturn(:)
    net_growth(:)           = net_growth(:) - lctlib_frac_sapwood_branch * veg_pool%sap_wood%carbon(:) * fturn(:)

    !> 4.2 heart wood formation, including nutrient recycling
    !! 
    IF(lctlib_growthform == itree)THEN
        fturn(:)              = (1.0_wp - lctlib_frac_sapwood_branch) / lctlib_tau_sap_wood / one_day / one_year * dtime
        veg_growth%heart_wood = veg_growth%heart_wood + veg_pool%sap_wood * fturn(:)
        veg_growth%sap_wood   = veg_growth%sap_wood   - veg_pool%sap_wood * fturn(:) 
        ! Nutrient recycling
        nitrogen_recycling_fract(:)   = calc_nitrogen_recycling_fract("wood")
        phosphorus_recycling_fract(:) = calc_phosphorus_recycling_fract("wood")
        recycling_heart_wood_n(:)     = recycling_heart_wood_n(:)   + & 
                                          veg_pool%sap_wood%nitrogen(:)   * fturn(:) * nitrogen_recycling_fract(:)
        recycling_heart_wood_n15(:)   = recycling_heart_wood_n15(:) + &
                                          veg_pool%sap_wood%nitrogen15(:) * fturn(:) * nitrogen_recycling_fract(:)
        recycling_heart_wood_p(:)     = recycling_heart_wood_p(:)   + & 
                                          veg_pool%sap_wood%phosphorus(:) * fturn(:) * phosphorus_recycling_fract(:)
    ENDIF 


    !> 5.0 fruit turnover to seed_bed
    !! 
    fturn(:)            = 1.0_wp / lctlib_tau_fruit / one_day / one_year * dtime
    veg_growth%seed_bed = veg_growth%seed_bed + veg_pool%fruit * fturn(:)
    veg_growth%fruit    = veg_growth%fruit    - veg_pool%fruit * fturn(:) 
    net_growth(:)       = net_growth(:)       - veg_pool%fruit%carbon(:) * fturn(:)

    !> 6.0 seed-bed turnover to litter 
    !! 
    fturn(:)                = 1.0_wp / lctlib_tau_seed_litter / one_day / one_year * dtime
    veg_litterfall%seed_bed = veg_litterfall%seed_bed + veg_pool%seed_bed * fturn(:)

  END SUBROUTINE calc_turnover


  !-----------------------------------------------------------------------------------------------------
  !> calculate fractional N recycling from leaves 
  !!
  !! should be made dynamically dependended on N stress later 
  !!
  !! @todo veg_turnover:calc_nitrogen_recycling_fract implement functionality other than "set to default"
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_nitrogen_recycling_fract(pool_name) RESULT (recycling_fract)

    USE mo_kind,                  ONLY: wp
    USE mo_veg_constants,         ONLY: resorp_fract_leaf, resorp_fract_wood

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(4), INTENT(in) :: pool_name
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                 :: recycling_fract 
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

    USE mo_kind,                  ONLY: wp
    USE mo_veg_constants,         ONLY: resorp_fract_leaf,resorp_fract_wood

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(4), INTENT(in) :: pool_name
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                 :: recycling_fract 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_phosphorus_recycling_fract'

    SELECT CASE(TRIM(pool_name))
      CASE("leaf")
        recycling_fract = resorp_fract_leaf 
      CASE("wood")
        recycling_fract = resorp_fract_wood
    ENDSELECT

  END FUNCTION calc_phosphorus_recycling_fract

END MODULE mo_veg_turnover
