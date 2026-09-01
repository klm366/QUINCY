!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### routine for updating fire fluxes
!>
MODULE mo_q_dfire_sf_feedback_class
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish

  USE mo_q_dfire_sf_constants, ONLY: t_sf_fuel_dead_fuel_woody_ix, t_sf_fuel_dead_fuel_polymeric_ix, & !< dead fuel pools indices (from sb_pool_mt)
                                      t_sf_fuel_dead_fuel_dom_assoc_ix, t_sf_fuel_dead_fuel_dom_ix, &
                                      t_sf_fuel_dead_fuel_res_ix, t_sf_fuel_dead_fuel_res_assoc_ix, &
                                      t_sf_fuel_dead_fuel_soluable_ix

  USE mo_q_dfire_sf_constants, ONLY: t_sf_fuel_live_fuel_leaf_ix, t_sf_fuel_live_fuel_sap_wood_ix, & !< live fuel pool indices
                                      t_sf_fuel_live_fuel_heart_wood_ix, t_sf_fuel_live_fuel_labile_ix, &
                                      t_sf_fuel_live_fuel_reserve_ix, t_sf_fuel_live_fuel_fruit_ix

  USE mo_q_dfire_sf_constants, ONLY: t_sf_fuel_dead_fuel_grass_fruit_ix, t_sf_fuel_dead_fuel_grass_labile_ix, & !< grass dead fuel pools indices (from veg_pool_mt)
                                      t_sf_fuel_dead_fuel_grass_leaf_ix, t_sf_fuel_dead_fuel_grass_reserve_ix, &
                                      t_sf_fuel_dead_fuel_grass_sap_ix

  USE mo_lnd_bgcm_idx,          ONLY: ix_woody_litter, ix_polymeric_litter, ix_soluable_litter, &
                                      ix_residue_assoc, ix_residue, &
                                      ix_dom_assoc, ix_dom, &
                                      ix_sap_wood, ix_heart_wood, ix_fruit, ix_labile, &
                                      ix_reserve, ix_leaf

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_sf_feedback_dead, calc_sf_feedback_live

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_feedback'

CONTAINS


  ! ======================================================================================================= !
  !>
  !> Apply spitfire effects into the model for dead fuel consumption on the top soil layer
  !>
  !> Woody (polymeric) litter pool (from sb_pool)
  !> -Emitted: add to sb_fire_emissions_mt (Removed)
  !> -Residue: add  to  sb_fire_litter_ash_mt -polymeric- (added)
  !>
  !> Live grass pools (from veg_pool):
  !> -Emitted: add to veg_fire_emissions_mt (substracted)
  !> -Residue: Add to veg_litterfall_mt (transfered from veg_pool into sb_pool)
  !>
  SUBROUTINE calc_sf_feedback_dead(dix, lix, &
      emissions, residue, &
      sb_fire_emissions_mt, sb_fire_litter_ash_mt, &
      veg_fire_emissions_mt, veg_litterfall_mt, &
      mask)

    INTEGER, INTENT(in)                   :: dix(:)                       !< Dead fuel elements indices (flag_sb_interactive_n, flag_sb_interactive_p)
    INTEGER, INTENT(in)                   :: lix(:)                       !< Live fuel elements indices (flag_veg_interactive_n, flag_veg_interactive_p)
    REAL(wp), INTENT(in)                  :: emissions(:,:,:)             !< Dead emissions (#nc, #dead fuel pools, #elems) [molX/m3 and molX/m2]
    REAL(wp), INTENT(in)                  :: residue(:,:,:)               !< Residue from dead fuel consumption (#nc, #dead fuel pools, #elems) [molX/m3 and molX/m2]
    REAL(wp), INTENT(inout)               :: veg_litterfall_mt(:,:,:)     !< bccm veg_pool (#pool, #element, #nc) [mol/m2]
    REAL(wp), INTENT(inout)               :: sb_fire_emissions_mt(:,:,:)  !< bgcm sb_pool (#pool, #element, #nc) [mol/m3]
    REAL(wp), INTENT(inout)               :: veg_fire_emissions_mt(:,:,:) !< bgcm veg_pool (#pool, #element, #nc) [mol/m2]
    REAL(wp), INTENT(inout)               :: sb_fire_litter_ash_mt(:,:,:) !< bgcm sb_pool (#pools, #element, #nc) [mol/m3]
    LOGICAL, INTENT(in)                   :: mask(:)                      !< True = Gridcells to compute

    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':calc_sf_feedback_dead'
    INTEGER :: inc
    ! ----------------------------------------------------------------------------------------------------- !
    !< Emitted FC is removed through sb_fire_emissions_mt
    DO inc=1, SIZE(sb_fire_emissions_mt, DIM=3) ! nc dim=3
      IF (mask(inc)) THEN
        sb_fire_emissions_mt(ix_woody_litter,dix,inc) = sb_fire_emissions_mt(ix_woody_litter,dix,inc) + &
                                                            emissions(inc,t_sf_fuel_dead_fuel_woody_ix,dix)
        sb_fire_emissions_mt(ix_polymeric_litter,dix,inc) = sb_fire_emissions_mt(ix_polymeric_litter,dix,inc) + &
                                                            emissions(inc,t_sf_fuel_dead_fuel_polymeric_ix,dix)
        sb_fire_emissions_mt(ix_soluable_litter,dix,inc) = sb_fire_emissions_mt(ix_soluable_litter,dix,inc) + &
                                                            emissions(inc,t_sf_fuel_dead_fuel_soluable_ix,dix)
        sb_fire_emissions_mt(ix_dom_assoc,dix,inc) = sb_fire_emissions_mt(ix_dom_assoc,dix,inc) + &
                                                            emissions(inc,t_sf_fuel_dead_fuel_dom_assoc_ix,dix)
        sb_fire_emissions_mt(ix_dom,dix,inc) = sb_fire_emissions_mt(ix_dom,dix,inc) + &
                                                            emissions(inc,t_sf_fuel_dead_fuel_dom_ix,dix)
        sb_fire_emissions_mt(ix_residue,dix,inc) = sb_fire_emissions_mt(ix_residue,dix,inc) + &
                                                            emissions(inc,t_sf_fuel_dead_fuel_res_ix,dix)
        sb_fire_emissions_mt(ix_residue_assoc,dix,inc) = sb_fire_emissions_mt(ix_residue_assoc,dix,inc) + &
                                                            emissions(inc,t_sf_fuel_dead_fuel_res_assoc_ix,dix)
      END IF
    END DO

    !> Emitted FC - livegrass (counts as dead fuel)
    DO inc=1, SIZE(veg_fire_emissions_mt, DIM=3) ! nc dim=3
      IF (mask(inc)) THEN
        veg_fire_emissions_mt(ix_fruit,lix,inc) = veg_fire_emissions_mt(ix_fruit,lix,inc) + &
            emissions(inc,t_sf_fuel_dead_fuel_grass_fruit_ix,lix)
        veg_fire_emissions_mt(ix_reserve,lix,inc) = veg_fire_emissions_mt(ix_reserve,lix,inc) + &
            emissions(inc,t_sf_fuel_dead_fuel_grass_reserve_ix,lix)
        veg_fire_emissions_mt(ix_labile,lix,inc) = veg_fire_emissions_mt(ix_labile,lix,inc) + &
            emissions(inc,t_sf_fuel_dead_fuel_grass_labile_ix,lix)
        veg_fire_emissions_mt(ix_sap_wood,lix,inc) = veg_fire_emissions_mt(ix_sap_wood,lix,inc) + &
            emissions(inc,t_sf_fuel_dead_fuel_grass_sap_ix,lix)
        veg_fire_emissions_mt(ix_leaf,lix,inc) = veg_fire_emissions_mt(ix_leaf,lix,inc) + &
            emissions(inc,t_sf_fuel_dead_fuel_grass_leaf_ix,lix)
      END IF
    END DO

    !< Residue is transfered using sb_fire_litter_ash_mt (woody litter to polymeric)
    DO inc=1, SIZE(sb_fire_litter_ash_mt, DIM=3) !< nc dim = 3
      IF (mask(inc)) THEN
        sb_fire_litter_ash_mt(ix_woody_litter,dix,inc) = &
                                          sb_fire_litter_ash_mt(ix_woody_litter,dix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_woody_ix,dix)
        sb_fire_litter_ash_mt(ix_polymeric_litter,dix,inc) = &
                                          sb_fire_litter_ash_mt(ix_polymeric_litter,dix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_polymeric_ix,dix)
        sb_fire_litter_ash_mt(ix_soluable_litter,dix,inc) = &
                                          sb_fire_litter_ash_mt(ix_soluable_litter,dix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_soluable_ix,dix)
        sb_fire_litter_ash_mt(ix_dom_assoc,dix,inc) = &
                                          sb_fire_litter_ash_mt(ix_dom_assoc,dix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_dom_assoc_ix,dix)
        sb_fire_litter_ash_mt(ix_dom,dix,inc) = &
                                          sb_fire_litter_ash_mt(ix_dom,dix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_dom_ix,dix)
        sb_fire_litter_ash_mt(ix_residue_assoc,dix,inc) = &
                                          sb_fire_litter_ash_mt(ix_residue_assoc,dix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_res_assoc_ix,dix)
        sb_fire_litter_ash_mt(ix_residue,dix,inc) = &
                                          sb_fire_litter_ash_mt(ix_residue,dix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_res_ix,dix)
        !< livegrass into veg_litterfall (no units conversion required)
        veg_litterfall_mt(ix_fruit,lix,inc) = veg_litterfall_mt(ix_fruit,lix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_grass_fruit_ix,lix)
        veg_litterfall_mt(ix_leaf,lix,inc) = veg_litterfall_mt(ix_leaf,lix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_grass_leaf_ix,lix)
        veg_litterfall_mt(ix_labile,lix,inc) = veg_litterfall_mt(ix_labile,lix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_grass_labile_ix,lix)
        veg_litterfall_mt(ix_sap_wood,lix,inc) = veg_litterfall_mt(ix_sap_wood,lix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_grass_sap_ix,lix)
        veg_litterfall_mt(ix_reserve,lix,inc) = veg_litterfall_mt(ix_reserve,lix,inc) + &
                                          residue(inc,t_sf_fuel_dead_fuel_grass_reserve_ix,lix)
      ENDIF
    END DO

  END SUBROUTINE calc_sf_feedback_dead


  ! ======================================================================================================= !
  !>
  !> Apply spitfire effects into the model for live fuel consumption
  !>
  !> Emitted: add to veg_fire_emissions (remove)
  !> Residue: add to veg_litterfall (transfer to polymeric litter-sb_pool-)
  !>
  SUBROUTINE calc_sf_feedback_live(lix, emissions, residue, veg_fire_emissions_mt, veg_litterfall_mt, mask)

    INTEGER, INTENT(in)                   :: lix(:)                       !< Live fuel element indices (flag_veg_interactive_n, flag_veg_interactive_p)
    REAL(wp), INTENT(in)                  :: emissions(:,:,:)             !< Live emissions (#nc,#n live pools,#elems) [mol/m2]
    REAL(wp), INTENT(in)                  :: residue(:,:,:)               !< Residue from live emissions (#nc,#n live pools,#elems) [mol/m2]
    REAL(wp), INTENT(inout)               :: veg_fire_emissions_mt(:,:,:) !< bgcm veg_pool (#pool, #element, #nc) [mol/m2]
    REAL(wp), INTENT(inout)               :: veg_litterfall_mt(:,:,:)     !< bgcm veg_pool (#pool, #element, #nc) [mol/m2]
    LOGICAL, INTENT(in)                   :: mask(:)                      !< True = compute gridcell(s)

    INTEGER                               :: inc
    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':calc_sf_feedback_live'
    ! ----------------------------------------------------------------------------------------------------- !
    !< Emitted
    DO inc=1, SIZE(emissions, DIM=1) !< nc dim=1
      IF (mask(inc)) THEN
        veg_fire_emissions_mt(ix_fruit,lix,inc) = veg_fire_emissions_mt(ix_fruit,lix,inc) + &
                                                emissions(inc,t_sf_fuel_live_fuel_fruit_ix,lix)
        veg_fire_emissions_mt(ix_leaf,lix,inc) = veg_fire_emissions_mt(ix_leaf,lix,inc) + &
                                                emissions(inc,t_sf_fuel_live_fuel_leaf_ix,lix)
        veg_fire_emissions_mt(ix_labile,lix,inc) = veg_fire_emissions_mt(ix_labile,lix,inc) + &
                                                emissions(inc,t_sf_fuel_live_fuel_labile_ix,lix)
        veg_fire_emissions_mt(ix_reserve,lix,inc) = veg_fire_emissions_mt(ix_reserve,lix,inc) + &
                                                emissions(inc,t_sf_fuel_live_fuel_reserve_ix,lix)
        veg_fire_emissions_mt(ix_sap_wood,lix,inc) = veg_fire_emissions_mt(ix_sap_wood,lix,inc) + &
                                                emissions(inc,t_sf_fuel_live_fuel_sap_wood_ix,lix)
        veg_fire_emissions_mt(ix_heart_wood,lix,inc) = veg_fire_emissions_mt(ix_heart_wood,lix,inc) + &
                                                emissions(inc,t_sf_fuel_live_fuel_heart_wood_ix,lix)
      END IF
    END DO

    !< Residue (add to veg_littefall_mt)
    DO inc=1, SIZE(residue, DIM=1) ! 1=nc dim
      IF (mask(inc)) THEN
        veg_litterfall_mt(ix_fruit,lix,inc) = veg_litterfall_mt(ix_fruit,lix,inc) + &
                                                      residue(inc,t_sf_fuel_live_fuel_fruit_ix,lix)
        veg_litterfall_mt(ix_leaf,lix,inc) = veg_litterfall_mt(ix_leaf,lix,inc) + &
                                                      residue(inc,t_sf_fuel_live_fuel_leaf_ix,lix)
        veg_litterfall_mt(ix_labile,lix,inc) = veg_litterfall_mt(ix_labile,lix,inc) + &
                                                      residue(inc,t_sf_fuel_live_fuel_labile_ix,lix)
        veg_litterfall_mt(ix_reserve,lix,inc) = veg_litterfall_mt(ix_reserve,lix,inc) + &
                                                      residue(inc,t_sf_fuel_live_fuel_reserve_ix,lix)
        veg_litterfall_mt(ix_sap_wood,lix,inc) = veg_litterfall_mt(ix_sap_wood,lix,inc) + &
                                                      residue(inc,t_sf_fuel_live_fuel_sap_wood_ix,lix)
        veg_litterfall_mt(ix_heart_wood,lix,inc) = veg_litterfall_mt(ix_heart_wood,lix,inc) + &
                                                      residue(inc,t_sf_fuel_live_fuel_heart_wood_ix,lix)
      END IF
    END DO

  END SUBROUTINE calc_sf_feedback_live
#endif

END MODULE mo_q_dfire_sf_feedback_class
