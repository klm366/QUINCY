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
!! @brief 
!!  Contains the science routines for the litter transport and partitioning for the JSM and SIMPLE soil models
!!
!! @author
!!  Sönke Zaehle 
!!
!! @par Revision History
!!  First version by Jan Engel (20181001): copyied from sb/mo_sb_simple_model.f90
!!
MODULE mo_sb_litter_processes
#ifndef __NO_JSBACH__

#ifdef __QUINCY_STANDALONE__
#include "dsl_jsbach.inc"
#endif

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_litter_partitioning

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_litter_processes'

CONTAINS

  ! ======================================================================================================= !
  !>calculate the partitioning of litter fall from vegetation into metabolic, structural, 
  !> and woody litter given: 
  !>
  !>  a) tissue specifc preferences (labile + reserve = metabolic sapwood and heartwood = woody)
  !>
  !>  b) the tissue-specific lignin concentration and the actual tissue C:N ratio.
  !>
  !>  All above-ground litter enters the first layer, only fine and coarse roots are distributed
  !>  according to their root profile 
  !>
  SUBROUTINE calc_litter_partitioning(nc, &                         ! in
                                      nsoil_sb, &
                                      lctlib_sla, &
                                      lctlib_growthform, &
                                      sb_model_scheme, &
                                      soil_depth_sl, &
                                      root_fraction_act_sl, &
                                      veg_litterfall, &
                                      sb_formation )                ! inout

    USE mo_veg_constants,           ONLY: itree, igrass
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_sb_constants,            ONLY: fc_soluable_max,k_fc_soluable,k_fn_soluable,k_fp_soluable, &
                                          lc_leaf_max,lc_leaf2sla,lc_fine_root,lc_coarse_root,lc_fruit,lc_seed_bed 
#ifdef __QUINCY_STANDALONE__
    dsl4slm_Use_memory_a_sb
    dsl4slm_Use_memory_a_veg
#else
    dsl4jsb_Use_pools(SB_)  t_sb_components_pool
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
#endif
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                      !< dimensions
                                                                       nsoil_sb
    REAL(wp),                                         INTENT(in)    :: lctlib_sla                 !< lctlib parameter
    INTEGER,                                          INTENT(in)    :: lctlib_growthform          !< lctlib parameter
    CHARACTER(len=*),                                 INTENT(in)    :: sb_model_scheme            !< soil model/module name
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl              !< depth per soil layer
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: root_fraction_act_sl       !< root fraction per soil layer
#ifdef __QUINCY_STANDALONE__
    TYPE(t_veg_pool_a),                               INTENT(in)    :: veg_litterfall             !< vegetation litterfall pool
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_formation               !< soil biogeochemistry formation pool
#else
    TYPE(t_veg_tissues_pool),                         INTENT(in)    :: veg_litterfall             !< vegetation litterfall pool
    TYPE(t_sb_components_pool),                       INTENT(inout) :: sb_formation               !< soil biogeochemistry formation pool
#endif
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)         :: fc_soluable_leaf, fc_soluable_fine_root, fc_soluable_coarse_root, fc_soluable_fruit, & 
                                       fc_soluable_seed_bed, fc_soluable_labile, fc_soluable_reserve 
    REAL(wp), DIMENSION(nc)         :: fn_soluable_leaf, fn_soluable_fine_root, fn_soluable_coarse_root, fn_soluable_fruit, & 
                                       fn_soluable_seed_bed, fn_soluable_labile, fn_soluable_reserve 
    REAL(wp), DIMENSION(nc)         :: fp_soluable_leaf, fp_soluable_fine_root, fp_soluable_coarse_root, fp_soluable_fruit, & 
                                       fp_soluable_seed_bed, fp_soluable_labile, fp_soluable_reserve 
    INTEGER                         :: ic                         ! loop over grid points
    INTEGER                         :: isoil                      ! loop over soil layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_partitioning'

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !>0.9 init local variables
    !>
    fc_soluable_leaf(:)         = 0.0_wp
    fn_soluable_leaf(:)         = 0.0_wp
    fp_soluable_leaf(:)         = 0.0_wp
    fc_soluable_fine_root(:)    = 0.0_wp 
    fn_soluable_fine_root(:)    = 0.0_wp 
    fp_soluable_fine_root(:)    = 0.0_wp 
    fc_soluable_coarse_root(:)  = 0.0_wp 
    fn_soluable_coarse_root(:)  = 0.0_wp 
    fp_soluable_coarse_root(:)  = 0.0_wp 
    fc_soluable_fruit(:)        = 0.0_wp
    fn_soluable_fruit(:)        = 0.0_wp
    fp_soluable_fruit(:)        = 0.0_wp
    fc_soluable_seed_bed(:)     = 0.0_wp
    fn_soluable_seed_bed(:)     = 0.0_wp
    fp_soluable_seed_bed(:)     = 0.0_wp
    fc_soluable_labile(:)       = 1.0_wp 
    fn_soluable_labile(:)       = 1.0_wp 
    fp_soluable_labile(:)       = 1.0_wp 
    fc_soluable_reserve(:)      = 1.0_wp
    fn_soluable_reserve(:)      = 1.0_wp
    fp_soluable_reserve(:)      = 1.0_wp

    !>1.0 determine allocation of fresh litter fall to litter pools
    !>
    WHERE(veg_litterfall%leaf%nitrogen(:) > eps8)
       fc_soluable_leaf(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * (lc_leaf_max + lc_leaf2sla * lctlib_sla) * & 
                             veg_litterfall%leaf%carbon(:) / veg_litterfall%leaf%nitrogen(:))
       WHERE(fc_soluable_leaf(:) > eps8)
          fn_soluable_leaf(:) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(:))/(k_fn_soluable*fc_soluable_leaf(:)))
          fp_soluable_leaf(:) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(:))/(k_fp_soluable*fc_soluable_leaf(:)))
       END WHERE
    END WHERE

    WHERE(veg_litterfall%fine_root%nitrogen(:) > eps8)
       fc_soluable_fine_root(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_fine_root * & 
                                  veg_litterfall%fine_root%carbon(:) / veg_litterfall%fine_root%nitrogen(:))
       WHERE(fc_soluable_fine_root(:) > eps8)
          fn_soluable_fine_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fine_root(:))/(k_fn_soluable*fc_soluable_fine_root(:)))
          fp_soluable_fine_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fine_root(:))/(k_fp_soluable*fc_soluable_fine_root(:)))
       END WHERE
    END WHERE

    WHERE(veg_litterfall%coarse_root%nitrogen(:) > eps8)
       fc_soluable_coarse_root(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_coarse_root * & 
                                    veg_litterfall%coarse_root%carbon(:) / veg_litterfall%coarse_root%nitrogen(:))
       WHERE(fc_soluable_coarse_root(:) > eps8)
          fn_soluable_coarse_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(:)) / &
                                       (k_fn_soluable*fc_soluable_coarse_root(:)))
          fp_soluable_coarse_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(:)) / &
                                       (k_fp_soluable*fc_soluable_coarse_root(:)))
        END WHERE
    END WHERE

    WHERE(veg_litterfall%fruit%nitrogen(:) > eps8)
       fc_soluable_fruit(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_fruit * & 
                              veg_litterfall%fruit%carbon(:) / veg_litterfall%fruit%nitrogen(:))
       WHERE(fc_soluable_fruit(:) > eps8)
          fn_soluable_fruit(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(:))/(k_fn_soluable*fc_soluable_fruit(:)))
          fp_soluable_fruit(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(:))/(k_fp_soluable*fc_soluable_fruit(:)))
       END WHERE
    END WHERE

    WHERE(veg_litterfall%seed_bed%nitrogen(:) > eps8)
       fc_soluable_seed_bed(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_seed_bed * & 
                                 veg_litterfall%seed_bed%carbon(:) / veg_litterfall%seed_bed%nitrogen(:))
       WHERE(fc_soluable_seed_bed(:) > eps8)
          fn_soluable_seed_bed(:) = 1._wp / (1._wp + (1._wp - fc_soluable_seed_bed(:))/(k_fn_soluable*fc_soluable_seed_bed(:)))
          fp_soluable_seed_bed(:) = 1._wp / (1._wp + (1._wp - fc_soluable_seed_bed(:))/(k_fp_soluable*fc_soluable_seed_bed(:)))
       END WHERE
    END WHERE


    ! loop over grid points and check for soil layer width > eps8
    DO ic = 1,nc
      IF (soil_depth_sl(ic,1) > eps8) THEN
        !>2.0 Update surface metabolic and soluable litter pools with litter fall
        !>
        sb_formation%soluable_litter%carbon(ic,1) = sb_formation%soluable_litter%carbon(ic,1) + & 
          & (fc_soluable_leaf(ic)     * veg_litterfall%leaf%carbon(ic)    + & 
          &  fc_soluable_fruit(ic)    * veg_litterfall%fruit%carbon(ic)   + & 
          &  fc_soluable_labile(ic)   * veg_litterfall%labile%carbon(ic)  + & 
          &  fc_soluable_reserve(ic)  * veg_litterfall%reserve%carbon(ic) + & 
          &  fc_soluable_seed_bed(ic) * veg_litterfall%seed_bed%carbon(ic) ) / soil_depth_sl(ic,1)
        sb_formation%soluable_litter%nitrogen(ic,1) = sb_formation%soluable_litter%nitrogen(ic,1) + & 
          & (fn_soluable_leaf(ic)     * veg_litterfall%leaf%nitrogen(ic)    + & 
          &  fn_soluable_fruit(ic)    * veg_litterfall%fruit%nitrogen(ic)   + & 
          &  fn_soluable_labile(ic)   * veg_litterfall%labile%nitrogen(ic)  + & 
          &  fn_soluable_reserve(ic)  * veg_litterfall%reserve%nitrogen(ic) + & 
          &  fn_soluable_seed_bed(ic) * veg_litterfall%seed_bed%nitrogen(ic) ) / soil_depth_sl(ic,1)
        sb_formation%soluable_litter%phosphorus(ic,1) = sb_formation%soluable_litter%phosphorus(ic,1) + & 
          & (fp_soluable_leaf(ic)     * veg_litterfall%leaf%phosphorus(ic)    + & 
          &  fp_soluable_fruit(ic)    * veg_litterfall%fruit%phosphorus(ic)   + & 
          &  fp_soluable_labile(ic)   * veg_litterfall%labile%phosphorus(ic)  + & 
          &  fp_soluable_reserve(ic)  * veg_litterfall%reserve%phosphorus(ic) + & 
          &  fp_soluable_seed_bed(ic) * veg_litterfall%seed_bed%phosphorus(ic) ) / soil_depth_sl(ic,1)
        sb_formation%soluable_litter%carbon13(ic,1) = sb_formation%soluable_litter%carbon13(ic,1) + & 
          & (fc_soluable_leaf(ic)     * veg_litterfall%leaf%carbon13(ic)    + & 
          &  fc_soluable_fruit(ic)    * veg_litterfall%fruit%carbon13(ic)   + & 
          &  fc_soluable_labile(ic)   * veg_litterfall%labile%carbon13(ic)  + & 
          &  fc_soluable_reserve(ic)  * veg_litterfall%reserve%carbon13(ic) + & 
          &  fc_soluable_seed_bed(ic) * veg_litterfall%seed_bed%carbon13(ic) ) / soil_depth_sl(ic,1)
        sb_formation%soluable_litter%carbon14(ic,1) = sb_formation%soluable_litter%carbon14(ic,1) + & 
          & (fc_soluable_leaf(ic)     * veg_litterfall%leaf%carbon14(ic)    + & 
          &  fc_soluable_fruit(ic)    * veg_litterfall%fruit%carbon14(ic)   + & 
          &  fc_soluable_labile(ic)   * veg_litterfall%labile%carbon14(ic)  + & 
          &  fc_soluable_reserve(ic)  * veg_litterfall%reserve%carbon14(ic) + & 
          &  fc_soluable_seed_bed(ic) * veg_litterfall%seed_bed%carbon14(ic) ) / soil_depth_sl(ic,1)
        sb_formation%soluable_litter%nitrogen15(ic,1) = sb_formation%soluable_litter%nitrogen15(ic,1) + & 
          & (fn_soluable_leaf(ic)     * veg_litterfall%leaf%nitrogen15(ic)    + & 
          &  fn_soluable_fruit(ic)    * veg_litterfall%fruit%nitrogen15(ic)   + & 
          &  fn_soluable_labile(ic)   * veg_litterfall%labile%nitrogen15(ic)  + & 
          &  fn_soluable_reserve(ic)  * veg_litterfall%reserve%nitrogen15(ic) + & 
          &  fn_soluable_seed_bed(ic) * veg_litterfall%seed_bed%nitrogen15(ic) ) / soil_depth_sl(ic,1)
        
        sb_formation%polymeric_litter%carbon(ic,1) = sb_formation%polymeric_litter%carbon(ic,1) + & 
          & ((1.0_wp - fc_soluable_leaf(ic))     * veg_litterfall%leaf%carbon(ic)    + & 
          &  (1.0_wp - fc_soluable_fruit(ic))    * veg_litterfall%fruit%carbon(ic)   + & 
          &  (1.0_wp - fc_soluable_labile(ic))   * veg_litterfall%labile%carbon(ic)  + & 
          &  (1.0_wp - fc_soluable_reserve(ic))  * veg_litterfall%reserve%carbon(ic) + & 
          &  (1.0_wp - fc_soluable_seed_bed(ic)) * veg_litterfall%seed_bed%carbon(ic) ) / soil_depth_sl(ic,1)
        sb_formation%polymeric_litter%nitrogen(ic,1) = sb_formation%polymeric_litter%nitrogen(ic,1) + & 
          & ((1.0_wp - fn_soluable_leaf(ic))     * veg_litterfall%leaf%nitrogen(ic)    + & 
          &  (1.0_wp - fn_soluable_fruit(ic))    * veg_litterfall%fruit%nitrogen(ic)   + & 
          &  (1.0_wp - fn_soluable_labile(ic))   * veg_litterfall%labile%nitrogen(ic)  + & 
          &  (1.0_wp - fn_soluable_reserve(ic))  * veg_litterfall%reserve%nitrogen(ic) + & 
          &  (1.0_wp - fn_soluable_seed_bed(ic)) * veg_litterfall%seed_bed%nitrogen(ic) ) / soil_depth_sl(ic,1)
        sb_formation%polymeric_litter%phosphorus(ic,1) = sb_formation%polymeric_litter%phosphorus(ic,1) + & 
          & ((1.0_wp - fp_soluable_leaf(ic))     * veg_litterfall%leaf%phosphorus(ic)    + & 
          &  (1.0_wp - fp_soluable_fruit(ic))    * veg_litterfall%fruit%phosphorus(ic)   + & 
          &  (1.0_wp - fp_soluable_labile(ic))   * veg_litterfall%labile%phosphorus(ic)  + & 
          &  (1.0_wp - fp_soluable_reserve(ic))  * veg_litterfall%reserve%phosphorus(ic) + & 
          &  (1.0_wp - fp_soluable_seed_bed(ic)) * veg_litterfall%seed_bed%phosphorus(ic) ) / soil_depth_sl(ic,1)
        sb_formation%polymeric_litter%carbon13(ic,1) = sb_formation%polymeric_litter%carbon13(ic,1) + & 
          & ((1.0_wp - fc_soluable_leaf(ic))     * veg_litterfall%leaf%carbon13(ic)    + & 
          &  (1.0_wp - fc_soluable_fruit(ic))    * veg_litterfall%fruit%carbon13(ic)   + & 
          &  (1.0_wp - fc_soluable_labile(ic))   * veg_litterfall%labile%carbon13(ic)  + & 
          &  (1.0_wp - fc_soluable_reserve(ic))  * veg_litterfall%reserve%carbon13(ic) + & 
          &  (1.0_wp - fc_soluable_seed_bed(ic)) * veg_litterfall%seed_bed%carbon13(ic) ) / soil_depth_sl(ic,1)
        sb_formation%polymeric_litter%carbon14(ic,1) = sb_formation%polymeric_litter%carbon14(ic,1) + & 
          & ((1.0_wp - fc_soluable_leaf(ic))     * veg_litterfall%leaf%carbon14(ic)    + & 
          &  (1.0_wp - fc_soluable_fruit(ic))    * veg_litterfall%fruit%carbon14(ic)   + & 
          &  (1.0_wp - fc_soluable_labile(ic))   * veg_litterfall%labile%carbon14(ic)  + & 
          &  (1.0_wp - fc_soluable_reserve(ic))  * veg_litterfall%reserve%carbon14(ic) + & 
          &  (1.0_wp - fc_soluable_seed_bed(ic)) * veg_litterfall%seed_bed%carbon14(ic) ) / soil_depth_sl(ic,1)
        sb_formation%polymeric_litter%nitrogen15(ic,1) = sb_formation%polymeric_litter%nitrogen15(ic,1) + & 
          & ((1.0_wp - fn_soluable_leaf(ic))     * veg_litterfall%leaf%nitrogen15(ic)    + & 
          &  (1.0_wp - fn_soluable_fruit(ic))    * veg_litterfall%fruit%nitrogen15(ic)   + & 
          &  (1.0_wp - fn_soluable_labile(ic))   * veg_litterfall%labile%nitrogen15(ic)  + & 
          &  (1.0_wp - fn_soluable_reserve(ic))  * veg_litterfall%reserve%nitrogen15(ic) + & 
          &  (1.0_wp - fn_soluable_seed_bed(ic)) * veg_litterfall%seed_bed%nitrogen15(ic) ) / soil_depth_sl(ic,1)

        !>3.0 Add surface litter fall from sapwood and heartwood pools to surface woody litter pool for trees
        !>     for grasses, stems and halms (=sapwood) are not actually wood, and therefore paritioned to 
        !>     metabolic and structural litter as if their were fine roots
        !>
        SELECT CASE(lctlib_growthform)
        CASE (itree)
          sb_formation%woody_litter%carbon(ic,1) = sb_formation%woody_litter%carbon(ic,1) + & 
            & (veg_litterfall%sap_wood%carbon(ic) + veg_litterfall%heart_wood%carbon(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%woody_litter%nitrogen(ic,1) = sb_formation%woody_litter%nitrogen(ic,1) +  & 
            & (veg_litterfall%sap_wood%nitrogen(ic) + veg_litterfall%heart_wood%nitrogen(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%woody_litter%phosphorus(ic,1) = sb_formation%woody_litter%phosphorus(ic,1) + & 
            & (veg_litterfall%sap_wood%phosphorus(ic) + veg_litterfall%heart_wood%phosphorus(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%woody_litter%carbon13(ic,1) = sb_formation%woody_litter%carbon13(ic,1) + & 
            & (veg_litterfall%sap_wood%carbon13(ic) + veg_litterfall%heart_wood%carbon13(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%woody_litter%carbon14(ic,1) = sb_formation%woody_litter%carbon14(ic,1) + & 
            & (veg_litterfall%sap_wood%carbon14(ic) + veg_litterfall%heart_wood%carbon14(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%woody_litter%nitrogen15(ic,1) = sb_formation%woody_litter%nitrogen15(ic,1) + &
            & (veg_litterfall%sap_wood%nitrogen15(ic) + veg_litterfall%heart_wood%nitrogen15(ic)) / &
            & soil_depth_sl(ic,1)
        CASE (igrass)
          sb_formation%soluable_litter%carbon(ic,1) = sb_formation%soluable_litter%carbon(ic,1) + & 
            & fc_soluable_fine_root(ic) * (veg_litterfall%sap_wood%carbon(ic) + veg_litterfall%heart_wood%carbon(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%soluable_litter%nitrogen(ic,1) = sb_formation%soluable_litter%nitrogen(ic,1) + & 
            & fn_soluable_fine_root(ic) * (veg_litterfall%sap_wood%nitrogen(ic) + veg_litterfall%heart_wood%nitrogen(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%soluable_litter%phosphorus(ic,1) = sb_formation%soluable_litter%phosphorus(ic,1) + & 
            & fp_soluable_fine_root(ic) * (veg_litterfall%sap_wood%phosphorus(ic) + veg_litterfall%heart_wood%phosphorus(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%soluable_litter%carbon13(ic,1) = sb_formation%soluable_litter%carbon13(ic,1) + & 
            & fc_soluable_fine_root(ic) * (veg_litterfall%sap_wood%carbon13(ic) + veg_litterfall%heart_wood%carbon13(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%soluable_litter%carbon14(ic,1) = sb_formation%soluable_litter%carbon14(ic,1) + & 
            & fc_soluable_fine_root(ic) * (veg_litterfall%sap_wood%carbon14(ic) + veg_litterfall%heart_wood%carbon14(ic)) / & 
            & soil_depth_sl(ic,1)
          sb_formation%soluable_litter%nitrogen15(ic,1) = sb_formation%soluable_litter%nitrogen15(ic,1) + & 
            & fn_soluable_fine_root(ic) * (veg_litterfall%sap_wood%nitrogen15(ic) + veg_litterfall%heart_wood%nitrogen15(ic)) / &
            & soil_depth_sl(ic,1)
        
          sb_formation%polymeric_litter%carbon(ic,1) = sb_formation%polymeric_litter%carbon(ic,1) + & 
            & (1.0_wp - fc_soluable_fine_root(ic)) * &
            & (veg_litterfall%sap_wood%carbon(ic) + veg_litterfall%heart_wood%carbon(ic)) / soil_depth_sl(ic,1)
          sb_formation%polymeric_litter%nitrogen(ic,1) = sb_formation%polymeric_litter%nitrogen(ic,1) + & 
            & (1.0_wp - fn_soluable_fine_root(ic)) * &
            & (veg_litterfall%sap_wood%nitrogen(ic) + veg_litterfall%heart_wood%nitrogen(ic)) / soil_depth_sl(ic,1)
          sb_formation%polymeric_litter%phosphorus(ic,1) = sb_formation%polymeric_litter%phosphorus(ic,1) + & 
            & (1.0_wp - fp_soluable_fine_root(ic)) * &
            & (veg_litterfall%sap_wood%phosphorus(ic) + veg_litterfall%heart_wood%phosphorus(ic)) / soil_depth_sl(ic,1)
          sb_formation%polymeric_litter%carbon13(ic,1) = sb_formation%polymeric_litter%carbon13(ic,1) + & 
            & (1.0_wp - fc_soluable_fine_root(ic)) * &
            & (veg_litterfall%sap_wood%carbon13(ic) + veg_litterfall%heart_wood%carbon13(ic)) / soil_depth_sl(ic,1)
          sb_formation%polymeric_litter%carbon14(ic,1) = sb_formation%polymeric_litter%carbon14(ic,1) + & 
            & (1.0_wp - fc_soluable_fine_root(ic)) * &
            & (veg_litterfall%sap_wood%carbon14(ic) + veg_litterfall%heart_wood%carbon14(ic)) / soil_depth_sl(ic,1)
          sb_formation%polymeric_litter%nitrogen15(ic,1) = sb_formation%polymeric_litter%nitrogen15(ic,1) + & 
            & (1.0_wp - fn_soluable_fine_root(ic)) * &
            & (veg_litterfall%sap_wood%nitrogen15(ic) + veg_litterfall%heart_wood%nitrogen15(ic)) / soil_depth_sl(ic,1)
        END SELECT

        !>4.0 Distribute fine and coarse root litter to soil depth according to the root profile 
        !>
        DO isoil = 1,nsoil_sb
          IF (soil_depth_sl(ic,isoil) > eps8) THEN
            sb_formation%soluable_litter%carbon(ic,isoil)     = sb_formation%soluable_litter%carbon(ic,isoil) + & 
              & (fc_soluable_fine_root(ic)   * veg_litterfall%fine_root%carbon(ic)     + &
              &  fc_soluable_coarse_root(ic) * veg_litterfall%coarse_root%carbon(ic) ) * & 
              &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
            sb_formation%soluable_litter%nitrogen(ic,isoil)   = sb_formation%soluable_litter%nitrogen(ic,isoil) + & 
              & (fn_soluable_fine_root(ic)   * veg_litterfall%fine_root%nitrogen(ic)     + &
              &  fn_soluable_coarse_root(ic) * veg_litterfall%coarse_root%nitrogen(ic) ) * & 
              &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
            sb_formation%soluable_litter%phosphorus(ic,isoil) = sb_formation%soluable_litter%phosphorus(ic,isoil) + & 
              & (fp_soluable_fine_root(ic)   * veg_litterfall%fine_root%phosphorus(ic)     + &
              &  fp_soluable_coarse_root(ic) * veg_litterfall%coarse_root%phosphorus(ic) ) * & 
              &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
            sb_formation%soluable_litter%carbon13(ic,isoil)   = sb_formation%soluable_litter%carbon13(ic,isoil) + & 
              & (fc_soluable_fine_root(ic)   * veg_litterfall%fine_root%carbon13(ic)     + &
              &  fc_soluable_coarse_root(ic) * veg_litterfall%coarse_root%carbon13(ic) ) * & 
              &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
            sb_formation%soluable_litter%carbon14(ic,isoil)   = sb_formation%soluable_litter%carbon14(ic,isoil) + & 
              & (fc_soluable_fine_root(ic)   * veg_litterfall%fine_root%carbon14(ic)     + &
              &  fc_soluable_coarse_root(ic) * veg_litterfall%coarse_root%carbon14(ic) ) * & 
              &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
            sb_formation%soluable_litter%nitrogen15(ic,isoil) = sb_formation%soluable_litter%nitrogen15(ic,isoil) + & 
              & (fn_soluable_fine_root(ic)   * veg_litterfall%fine_root%nitrogen15(ic)     + &
              &  fn_soluable_coarse_root(ic) * veg_litterfall%coarse_root%nitrogen15(ic) ) * & 
              &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 

            IF (sb_model_scheme == "simple_1d") THEN
              sb_formation%polymeric_litter%carbon(ic,isoil)    = sb_formation%polymeric_litter%carbon(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_fine_root(ic))  * veg_litterfall%fine_root%carbon(ic)     + &
                &  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%carbon(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%nitrogen(ic,isoil)    = sb_formation%polymeric_litter%nitrogen(ic,isoil) + & 
                & ((1.0_wp - fn_soluable_fine_root(ic))  * veg_litterfall%fine_root%nitrogen(ic)     + &
                &  (1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%nitrogen(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%phosphorus(ic,isoil)    = sb_formation%polymeric_litter%phosphorus(ic,isoil) + & 
                & ((1.0_wp - fp_soluable_fine_root(ic))  * veg_litterfall%fine_root%phosphorus(ic)     + &
                &  (1.0_wp - fp_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%phosphorus(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%carbon13(ic,isoil)    = sb_formation%polymeric_litter%carbon13(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_fine_root(ic))  * veg_litterfall%fine_root%carbon13(ic)     + &
                &  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%carbon13(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%carbon14(ic,isoil)    = sb_formation%polymeric_litter%carbon14(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_fine_root(ic))  * veg_litterfall%fine_root%carbon14(ic)     + &
                &  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%carbon14(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%nitrogen15(ic,isoil)    = sb_formation%polymeric_litter%nitrogen15(ic,isoil) + & 
                & ((1.0_wp - fn_soluable_fine_root(ic))  * veg_litterfall%fine_root%nitrogen15(ic)     + &
                &  (1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%nitrogen15(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
            ELSE
              sb_formation%polymeric_litter%carbon(ic,isoil)    = sb_formation%polymeric_litter%carbon(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_fine_root(ic))  * veg_litterfall%fine_root%carbon(ic)) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%nitrogen(ic,isoil)    = sb_formation%polymeric_litter%nitrogen(ic,isoil) + & 
                & ((1.0_wp - fn_soluable_fine_root(ic))  * veg_litterfall%fine_root%nitrogen(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%phosphorus(ic,isoil)    = sb_formation%polymeric_litter%phosphorus(ic,isoil) + & 
                & ((1.0_wp - fp_soluable_fine_root(ic))  * veg_litterfall%fine_root%phosphorus(ic)  ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%carbon13(ic,isoil)    = sb_formation%polymeric_litter%carbon13(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_fine_root(ic))  * veg_litterfall%fine_root%carbon13(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%carbon14(ic,isoil)    = sb_formation%polymeric_litter%carbon14(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_fine_root(ic))  * veg_litterfall%fine_root%carbon14(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%polymeric_litter%nitrogen15(ic,isoil)    = sb_formation%polymeric_litter%nitrogen15(ic,isoil) + & 
                & ((1.0_wp - fn_soluable_fine_root(ic))  * veg_litterfall%fine_root%nitrogen15(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 

              sb_formation%woody_litter%carbon(ic,isoil)    = sb_formation%woody_litter%carbon(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%carbon(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%woody_litter%nitrogen(ic,isoil)    = sb_formation%woody_litter%nitrogen(ic,isoil) + & 
                & ((1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%nitrogen(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%woody_litter%phosphorus(ic,isoil)    = sb_formation%woody_litter%phosphorus(ic,isoil) + & 
                & ((1.0_wp - fp_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%phosphorus(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%woody_litter%carbon13(ic,isoil)    = sb_formation%woody_litter%carbon13(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%carbon13(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%woody_litter%carbon14(ic,isoil)    = sb_formation%woody_litter%carbon14(ic,isoil) + & 
                & ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%carbon14(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
              sb_formation%woody_litter%nitrogen15(ic,isoil)    = sb_formation%woody_litter%nitrogen15(ic,isoil) + & 
                & ((1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall%coarse_root%nitrogen15(ic) ) * & 
                &  root_fraction_act_sl(ic,isoil) / soil_depth_sl(ic,isoil) 
            END IF ! IF soil_depth_sl(ic,isoil) > eps8
          END IF   ! IF "simple_1d"
        END DO ! isoil = 1,nsoil_sb
      END IF   ! IF soil_depth_sl(ic,1) > eps8
    END DO     ! ic = 1,nc

  END SUBROUTINE calc_litter_partitioning

#endif
END MODULE mo_sb_litter_processes
