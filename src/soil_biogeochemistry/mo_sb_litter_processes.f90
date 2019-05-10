!>
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

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_litter_partitioning


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_litter_processes'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model and update_sb_jsm
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the partitioning of litter fall from vegetation into metabolic, structural, 
  !! and woody litter given a) tissue specifc preferences (labile + reserve = metabolic; 
  !! sapwood and heartwood = woody); b) the tissue-specific lignin concentration and the actual
  !! tissue C:N ratio.
  !! All above-ground litter enters the first layer, only fine and coarse roots are distributed
  !! according to their root profile 
  !!
  !-----------------------------------------------------------------------------------------------------
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
    dsl4slm_Use_memory_a_veg
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)    :: nc, &                      !< dimensions
                                                                       nsoil_sb
    REAL(wp),                                         INTENT(in)    :: lctlib_sla                 !< lctlib parameter
    INTEGER,                                          INTENT(in)    :: lctlib_growthform          !< lctlib parameter
    CHARACTER(len=*),                                 INTENT(in)    :: sb_model_scheme            !< soil model/module name
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl, &           !< depth per soil layer
                                                                       root_fraction_act_sl       !< root fraction per soil layer
    TYPE(t_veg_pool_a),                               INTENT(in)    :: veg_litterfall             !< vegetation litterfall pool
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_formation               !< soil biogeochemistry formation pool
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)         :: fc_soluable_leaf, fc_soluable_fine_root, fc_soluable_coarse_root, fc_soluable_fruit, & 
                                       fc_soluable_seed_bed, fc_soluable_labile, fc_soluable_reserve 
    REAL(wp), DIMENSION(nc)         :: fn_soluable_leaf, fn_soluable_fine_root, fn_soluable_coarse_root, fn_soluable_fruit, & 
                                       fn_soluable_seed_bed, fn_soluable_labile, fn_soluable_reserve 
    REAL(wp), DIMENSION(nc)         :: fp_soluable_leaf, fp_soluable_fine_root, fp_soluable_coarse_root, fp_soluable_fruit, & 
                                       fp_soluable_seed_bed, fp_soluable_labile, fp_soluable_reserve 
    INTEGER                         :: isoil                      !< loop over soil layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_partitioning'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 determine allocation of fresh litter fall to litter pools
    !!
    fc_soluable_leaf(:) = 0.0_wp
    fn_soluable_leaf(:) = 0.0_wp
    fp_soluable_leaf(:) = 0.0_wp
    WHERE(veg_litterfall%leaf%nitrogen(:) > eps8)
       fc_soluable_leaf(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * (lc_leaf_max + lc_leaf2sla * lctlib_sla) * & 
                             veg_litterfall%leaf%carbon(:) / veg_litterfall%leaf%nitrogen(:))
       WHERE(fc_soluable_leaf(:) > eps8)
          fn_soluable_leaf(:) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(:))/(k_fn_soluable*fc_soluable_leaf(:)))
          fp_soluable_leaf(:) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(:))/(k_fp_soluable*fc_soluable_leaf(:)))
       END WHERE
    END WHERE

    fc_soluable_fine_root(:) = 0.0_wp 
    fn_soluable_fine_root(:) = 0.0_wp 
    fp_soluable_fine_root(:) = 0.0_wp 
    WHERE(veg_litterfall%fine_root%nitrogen(:) > eps8)
       fc_soluable_fine_root(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_fine_root * & 
                                  veg_litterfall%fine_root%carbon(:) / veg_litterfall%fine_root%nitrogen(:))
       WHERE(fc_soluable_fine_root(:) > eps8)
          fn_soluable_fine_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fine_root(:))/(k_fn_soluable*fc_soluable_fine_root(:)))
          fp_soluable_fine_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fine_root(:))/(k_fp_soluable*fc_soluable_fine_root(:)))
       END WHERE
    END WHERE

    fc_soluable_coarse_root(:) = 0.0_wp 
    fn_soluable_coarse_root(:) = 0.0_wp 
    fp_soluable_coarse_root(:) = 0.0_wp 
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

    fc_soluable_fruit(:) = 0.0_wp
    fn_soluable_fruit(:) = 0.0_wp
    fp_soluable_fruit(:) = 0.0_wp
    WHERE(veg_litterfall%fruit%nitrogen(:) > eps8)
       fc_soluable_fruit(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_fruit * & 
                              veg_litterfall%fruit%carbon(:) / veg_litterfall%fruit%nitrogen(:))
       WHERE(fc_soluable_fruit(:) > eps8)
          fn_soluable_fruit(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(:))/(k_fn_soluable*fc_soluable_fruit(:)))
          fp_soluable_fruit(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(:))/(k_fp_soluable*fc_soluable_fruit(:)))
       END WHERE
    END WHERE

    fc_soluable_seed_bed(:) = 0.0_wp
    fn_soluable_seed_bed(:) = 0.0_wp
    fp_soluable_seed_bed(:) = 0.0_wp
    WHERE(veg_litterfall%seed_bed%nitrogen(:) > eps8)
       fc_soluable_seed_bed(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_seed_bed * & 
                                 veg_litterfall%seed_bed%carbon(:) / veg_litterfall%seed_bed%nitrogen(:))
       WHERE(fc_soluable_seed_bed(:) > eps8)
          fn_soluable_seed_bed(:) = 1._wp / (1._wp + (1._wp - fc_soluable_seed_bed(:))/(k_fn_soluable*fc_soluable_seed_bed(:)))
          fp_soluable_seed_bed(:) = 1._wp / (1._wp + (1._wp - fc_soluable_seed_bed(:))/(k_fp_soluable*fc_soluable_seed_bed(:)))
       END WHERE
    END WHERE

    fc_soluable_labile(:)  = 1.0_wp 
    fn_soluable_labile(:)  = 1.0_wp 
    fp_soluable_labile(:)  = 1.0_wp 
    fc_soluable_reserve(:) = 1.0_wp
    fn_soluable_reserve(:) = 1.0_wp
    fp_soluable_reserve(:) = 1.0_wp

    !> 2.0 Update surface metabolic and soluable litter pools with litter fall
    !!
    sb_formation%soluable_litter%carbon(:,1) = sb_formation%soluable_litter%carbon(:,1) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon(:)    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon(:)   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon(:)  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon(:) + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%nitrogen(:,1) = sb_formation%soluable_litter%nitrogen(:,1) + & 
       (fn_soluable_leaf(:)     * veg_litterfall%leaf%nitrogen(:)    + & 
        fn_soluable_fruit(:)    * veg_litterfall%fruit%nitrogen(:)   + & 
        fn_soluable_labile(:)   * veg_litterfall%labile%nitrogen(:)  + & 
        fn_soluable_reserve(:)  * veg_litterfall%reserve%nitrogen(:) + & 
        fn_soluable_seed_bed(:) * veg_litterfall%seed_bed%nitrogen(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%phosphorus(:,1) = sb_formation%soluable_litter%phosphorus(:,1) + & 
       (fp_soluable_leaf(:)     * veg_litterfall%leaf%phosphorus(:)    + & 
        fp_soluable_fruit(:)    * veg_litterfall%fruit%phosphorus(:)   + & 
        fp_soluable_labile(:)   * veg_litterfall%labile%phosphorus(:)  + & 
        fp_soluable_reserve(:)  * veg_litterfall%reserve%phosphorus(:) + & 
        fp_soluable_seed_bed(:) * veg_litterfall%seed_bed%phosphorus(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%carbon13(:,1) = sb_formation%soluable_litter%carbon13(:,1) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon13(:)    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon13(:)   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon13(:)  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon13(:) + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon13(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%carbon14(:,1) = sb_formation%soluable_litter%carbon14(:,1) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon14(:)    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon14(:)   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon14(:)  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon14(:) + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon14(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%nitrogen15(:,1) = sb_formation%soluable_litter%nitrogen15(:,1) + & 
       (fn_soluable_leaf(:)     * veg_litterfall%leaf%nitrogen15(:)    + & 
        fn_soluable_fruit(:)    * veg_litterfall%fruit%nitrogen15(:)   + & 
        fn_soluable_labile(:)   * veg_litterfall%labile%nitrogen15(:)  + & 
        fn_soluable_reserve(:)  * veg_litterfall%reserve%nitrogen15(:) + & 
        fn_soluable_seed_bed(:) * veg_litterfall%seed_bed%nitrogen15(:) ) / soil_depth_sl(:,1)

    sb_formation%polymeric_litter%carbon(:,1) = sb_formation%polymeric_litter%carbon(:,1) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon(:)    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon(:)   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon(:)  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon(:) + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%nitrogen(:,1) = sb_formation%polymeric_litter%nitrogen(:,1) + & 
       ((1.0_wp - fn_soluable_leaf(:))     * veg_litterfall%leaf%nitrogen(:)    + & 
        (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall%fruit%nitrogen(:)   + & 
        (1.0_wp - fn_soluable_labile(:))   * veg_litterfall%labile%nitrogen(:)  + & 
        (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall%reserve%nitrogen(:) + & 
        (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall%seed_bed%nitrogen(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%phosphorus(:,1) = sb_formation%polymeric_litter%phosphorus(:,1) + & 
       ((1.0_wp - fp_soluable_leaf(:))     * veg_litterfall%leaf%phosphorus(:)    + & 
        (1.0_wp - fp_soluable_fruit(:))    * veg_litterfall%fruit%phosphorus(:)   + & 
        (1.0_wp - fp_soluable_labile(:))   * veg_litterfall%labile%phosphorus(:)  + & 
        (1.0_wp - fp_soluable_reserve(:))  * veg_litterfall%reserve%phosphorus(:) + & 
        (1.0_wp - fp_soluable_seed_bed(:)) * veg_litterfall%seed_bed%phosphorus(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%carbon13(:,1) = sb_formation%polymeric_litter%carbon13(:,1) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon13(:)    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon13(:)   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon13(:)  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon13(:) + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon13(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%carbon14(:,1) = sb_formation%polymeric_litter%carbon14(:,1) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon14(:)    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon14(:)   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon14(:)  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon14(:) + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon14(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%nitrogen15(:,1) = sb_formation%polymeric_litter%nitrogen15(:,1) + & 
       ((1.0_wp - fn_soluable_leaf(:))     * veg_litterfall%leaf%nitrogen15(:)    + & 
        (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall%fruit%nitrogen15(:)   + & 
        (1.0_wp - fn_soluable_labile(:))   * veg_litterfall%labile%nitrogen15(:)  + & 
        (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall%reserve%nitrogen15(:) + & 
        (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall%seed_bed%nitrogen15(:) ) / soil_depth_sl(:,1)

    !> 3.0 Add surface litter fall from sapwood and heartwood pools to surface woody litter pool for trees
    !!     for grasses, stems and halms (=sapwood) are not actually wood, and therefore paritioned to 
    !!     metabolic and structural litter as if their were fine roots
    SELECT CASE(lctlib_growthform)
       CASE (itree)
          sb_formation%woody_litter%carbon(:,1) = sb_formation%woody_litter%carbon(:,1) + & 
             (veg_litterfall%sap_wood%carbon(:) + veg_litterfall%heart_wood%carbon(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%nitrogen(:,1) = sb_formation%woody_litter%nitrogen(:,1) +  & 
             (veg_litterfall%sap_wood%nitrogen(:) + veg_litterfall%heart_wood%nitrogen(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%phosphorus(:,1) = sb_formation%woody_litter%phosphorus(:,1) + & 
             (veg_litterfall%sap_wood%phosphorus(:) + veg_litterfall%heart_wood%phosphorus(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%carbon13(:,1) = sb_formation%woody_litter%carbon13(:,1) + & 
             (veg_litterfall%sap_wood%carbon13(:) + veg_litterfall%heart_wood%carbon13(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%carbon14(:,1) = sb_formation%woody_litter%carbon14(:,1) + & 
             (veg_litterfall%sap_wood%carbon14(:) + veg_litterfall%heart_wood%carbon14(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%nitrogen15(:,1) = sb_formation%woody_litter%nitrogen15(:,1) + &
             (veg_litterfall%sap_wood%nitrogen15(:) + veg_litterfall%heart_wood%nitrogen15(:)) / &
             soil_depth_sl(:,1)
       CASE (igrass)
          sb_formation%soluable_litter%carbon(:,1) = sb_formation%soluable_litter%carbon(:,1) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon(:) + veg_litterfall%heart_wood%carbon(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%nitrogen(:,1) = sb_formation%soluable_litter%nitrogen(:,1) + & 
             fn_soluable_fine_root(:) * (veg_litterfall%sap_wood%nitrogen(:) + veg_litterfall%heart_wood%nitrogen(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%phosphorus(:,1) = sb_formation%soluable_litter%phosphorus(:,1) + & 
             fp_soluable_fine_root(:) * (veg_litterfall%sap_wood%phosphorus(:) + veg_litterfall%heart_wood%phosphorus(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%carbon13(:,1) = sb_formation%soluable_litter%carbon13(:,1) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon13(:) + veg_litterfall%heart_wood%carbon13(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%carbon14(:,1) = sb_formation%soluable_litter%carbon14(:,1) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon14(:) + veg_litterfall%heart_wood%carbon14(:)) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%nitrogen15(:,1) = sb_formation%soluable_litter%nitrogen15(:,1) + & 
             fn_soluable_fine_root(:) * (veg_litterfall%sap_wood%nitrogen15(:) + veg_litterfall%heart_wood%nitrogen15(:)) / &
             soil_depth_sl(:,1)
      
          sb_formation%polymeric_litter%carbon(:,1) = sb_formation%polymeric_litter%carbon(:,1) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * &
             (veg_litterfall%sap_wood%carbon(:) + veg_litterfall%heart_wood%carbon(:)) / soil_depth_sl(:,1)
          sb_formation%polymeric_litter%nitrogen(:,1) = sb_formation%polymeric_litter%nitrogen(:,1) + & 
             (1.0_wp - fn_soluable_fine_root(:)) * &
             (veg_litterfall%sap_wood%nitrogen(:) + veg_litterfall%heart_wood%nitrogen(:)) / soil_depth_sl(:,1)
          sb_formation%polymeric_litter%phosphorus(:,1) = sb_formation%polymeric_litter%phosphorus(:,1) + & 
             (1.0_wp - fp_soluable_fine_root(:)) * &
             (veg_litterfall%sap_wood%phosphorus(:) + veg_litterfall%heart_wood%phosphorus(:)) / soil_depth_sl(:,1)
          sb_formation%polymeric_litter%carbon13(:,1) = sb_formation%polymeric_litter%carbon13(:,1) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * &
             (veg_litterfall%sap_wood%carbon13(:) + veg_litterfall%heart_wood%carbon13(:)) / soil_depth_sl(:,1)
          sb_formation%polymeric_litter%carbon14(:,1) = sb_formation%polymeric_litter%carbon14(:,1) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * &
             (veg_litterfall%sap_wood%carbon14(:) + veg_litterfall%heart_wood%carbon14(:)) / soil_depth_sl(:,1)
          sb_formation%polymeric_litter%nitrogen15(:,1) = sb_formation%polymeric_litter%nitrogen15(:,1) + & 
             (1.0_wp - fn_soluable_fine_root(:)) * &
             (veg_litterfall%sap_wood%nitrogen15(:) + veg_litterfall%heart_wood%nitrogen15(:)) / soil_depth_sl(:,1)

    END SELECT

    !> 4.0 Distribute fine and coarse root litter to soil depth according to the root profile 
    !!
    DO isoil = 1,nsoil_sb
       sb_formation%soluable_litter%carbon(:,isoil)     = sb_formation%soluable_litter%carbon(:,isoil) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon(:)     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%nitrogen(:,isoil)   = sb_formation%soluable_litter%nitrogen(:,isoil) + & 
          (fn_soluable_fine_root(:)   * veg_litterfall%fine_root%nitrogen(:)     + &
           fn_soluable_coarse_root(:) * veg_litterfall%coarse_root%nitrogen(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%phosphorus(:,isoil) = sb_formation%soluable_litter%phosphorus(:,isoil) + & 
          (fp_soluable_fine_root(:)   * veg_litterfall%fine_root%phosphorus(:)     + &
           fp_soluable_coarse_root(:) * veg_litterfall%coarse_root%phosphorus(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%carbon13(:,isoil)   = sb_formation%soluable_litter%carbon13(:,isoil) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon13(:)     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon13(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%carbon14(:,isoil)   = sb_formation%soluable_litter%carbon14(:,isoil) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon14(:)     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon14(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%nitrogen15(:,isoil) = sb_formation%soluable_litter%nitrogen15(:,isoil) + & 
          (fn_soluable_fine_root(:)   * veg_litterfall%fine_root%nitrogen15(:)     + &
           fn_soluable_coarse_root(:) * veg_litterfall%coarse_root%nitrogen15(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 

       IF (sb_model_scheme/="jsm") THEN
          sb_formation%polymeric_litter%carbon(:,isoil)    = sb_formation%polymeric_litter%carbon(:,isoil) + & 
             ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon(:)     + &
              (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%nitrogen(:,isoil)    = sb_formation%polymeric_litter%nitrogen(:,isoil) + & 
             ((1.0_wp - fn_soluable_fine_root(:))  * veg_litterfall%fine_root%nitrogen(:)     + &
              (1.0_wp - fn_soluable_coarse_root(:)) * veg_litterfall%coarse_root%nitrogen(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%phosphorus(:,isoil)    = sb_formation%polymeric_litter%phosphorus(:,isoil) + & 
             ((1.0_wp - fp_soluable_fine_root(:))  * veg_litterfall%fine_root%phosphorus(:)     + &
              (1.0_wp - fp_soluable_coarse_root(:)) * veg_litterfall%coarse_root%phosphorus(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%carbon13(:,isoil)    = sb_formation%polymeric_litter%carbon13(:,isoil) + & 
             ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon13(:)     + &
              (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon13(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%carbon14(:,isoil)    = sb_formation%polymeric_litter%carbon14(:,isoil) + & 
             ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon14(:)     + &
              (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon14(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%nitrogen15(:,isoil)    = sb_formation%polymeric_litter%nitrogen15(:,isoil) + & 
             ((1.0_wp - fn_soluable_fine_root(:))  * veg_litterfall%fine_root%nitrogen15(:)     + &
              (1.0_wp - fn_soluable_coarse_root(:)) * veg_litterfall%coarse_root%nitrogen15(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       ELSE
          sb_formation%polymeric_litter%carbon(:,isoil)    = sb_formation%polymeric_litter%carbon(:,isoil) + & 
             ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon(:)) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%nitrogen(:,isoil)    = sb_formation%polymeric_litter%nitrogen(:,isoil) + & 
             ((1.0_wp - fn_soluable_fine_root(:))  * veg_litterfall%fine_root%nitrogen(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%phosphorus(:,isoil)    = sb_formation%polymeric_litter%phosphorus(:,isoil) + & 
             ((1.0_wp - fp_soluable_fine_root(:))  * veg_litterfall%fine_root%phosphorus(:)  ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%carbon13(:,isoil)    = sb_formation%polymeric_litter%carbon13(:,isoil) + & 
             ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon13(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%carbon14(:,isoil)    = sb_formation%polymeric_litter%carbon14(:,isoil) + & 
             ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon14(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%polymeric_litter%nitrogen15(:,isoil)    = sb_formation%polymeric_litter%nitrogen15(:,isoil) + & 
             ((1.0_wp - fn_soluable_fine_root(:))  * veg_litterfall%fine_root%nitrogen15(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 

          sb_formation%woody_litter%carbon(:,isoil)    = sb_formation%woody_litter%carbon(:,isoil) + & 
             ((1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%woody_litter%nitrogen(:,isoil)    = sb_formation%woody_litter%nitrogen(:,isoil) + & 
             ((1.0_wp - fn_soluable_coarse_root(:)) * veg_litterfall%coarse_root%nitrogen(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%woody_litter%phosphorus(:,isoil)    = sb_formation%woody_litter%phosphorus(:,isoil) + & 
             ((1.0_wp - fp_soluable_coarse_root(:)) * veg_litterfall%coarse_root%phosphorus(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%woody_litter%carbon13(:,isoil)    = sb_formation%woody_litter%carbon13(:,isoil) + & 
             ((1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon13(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%woody_litter%carbon14(:,isoil)    = sb_formation%woody_litter%carbon14(:,isoil) + & 
             ((1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon14(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
          sb_formation%woody_litter%nitrogen15(:,isoil)    = sb_formation%woody_litter%nitrogen15(:,isoil) + & 
             ((1.0_wp - fn_soluable_coarse_root(:)) * veg_litterfall%coarse_root%nitrogen15(:) ) * & 
              root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       ENDIF
    END DO  

  END SUBROUTINE calc_litter_partitioning






END MODULE mo_sb_litter_processes
