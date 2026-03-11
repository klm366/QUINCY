!> QUINCY calculate fuel pools
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!> @description It holds and calculates the fuel from vegetation and litter.
!>
!>  Data is organized to account for:
!>    - Live fuel: from tree
!>        - Pool of origin
!>    - Dead fuel: from litter pool + 1 hr live grass, following Thonicke et al. (2010)
!>        - Pool of origin
!>    - Nutrients (C)
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_sf_fuel_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: finish, message_text

  USE mo_veg_constants,       ONLY: igrass, itree
  !< Top soil layers from sb_pool_mt used by spitfire, Total number of fuel classes
  USE mo_q_dfire_constants,   ONLY: df_sf_fc_full_nr
  !< Total number of spitfire dead fuel pools, polymeric fuel pool index within spitfire, same for woody litter
  USE mo_q_dfire_constants,   ONLY: sf_fuel_dead_fuel_nr, t_sf_fuel_dead_fuel_polymeric_ix, &
                                    t_sf_fuel_dead_fuel_woody_ix, t_sf_fuel_dead_fuel_dom_assoc_ix, &
                                    t_sf_fuel_dead_fuel_dom_ix, t_sf_fuel_dead_fuel_soluable_ix, &
                                    t_sf_fuel_dead_fuel_res_ix, t_sf_fuel_dead_fuel_res_assoc_ix
  !< Total number of live fuel pools within spitfire, live fuel pool indices from spitfire
  USE mo_q_dfire_constants,   ONLY: sf_fuel_live_fuel_nr, t_sf_fuel_live_fuel_sap_wood_ix, &
                                        t_sf_fuel_live_fuel_heart_wood_ix, t_sf_fuel_live_fuel_labile_ix, &
                                        t_sf_fuel_live_fuel_reserve_ix, t_sf_fuel_live_fuel_fruit_ix, &
                                        t_sf_fuel_live_fuel_leaf_ix
  !< Dead fuel pool indices from spitfire (live grass)
  USE mo_q_dfire_constants,   ONLY: t_sf_fuel_dead_fuel_grass_leaf_ix, t_sf_fuel_dead_fuel_grass_sap_ix, &
                                      t_sf_fuel_dead_fuel_grass_labile_ix, t_sf_fuel_dead_fuel_grass_reserve_ix, &
                                      t_sf_fuel_dead_fuel_grass_fruit_ix
  USE mo_q_dfire_constants,   ONLY: t_sf_litter_ratios
  !< Quincy sb_pool_mt pool indices
  USE mo_lnd_bgcm_idx,        ONLY: ix_woody_litter, ix_polymeric_litter, ix_leaf, ix_fruit, &
                                      ix_sap_wood, ix_heart_wood, ix_labile, ix_reserve, &
                                      ix_dom_assoc, ix_dom, ix_soluable_litter, ix_residue, ix_residue_assoc

  IMPLICIT NONE


  PUBLIC :: t_sf_fuel_class, sf_calc_has_any_fuel, sf_calc_aggr_fuels
  PUBLIC :: t_sf_fuel_dead_fuel_woody_ix, t_sf_fuel_dead_fuel_polymeric_ix
  PUBLIC :: create_fuel_class ! DEBUG only

  PRIVATE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_fuel_class'

  INTERFACE sf_calc_has_any_fuel
    MODULE PROCEDURE :: sf_calc_has_any_fuel_r2, sf_calc_has_any_fuel_r3
  END INTERFACE

  !< Spitfire fuel pool class
  TYPE :: t_sf_fuel_class
    INTEGER :: nc
    INTEGER, ALLOCATABLE :: growth_form(:)    !< itree, igrass ! todo: remove, not necessary
    TYPE(t_sf_litter_ratios) :: litter_ratios !< Litter ratio [0 to 1]
    REAL(wp), ALLOCATABLE :: cured_lvg_fract(:)!< Fraction of live grass considered dead fuel [0 to 1], Thonicke=1, Oberhagemann=0 to 1
    INTEGER :: ndeadpools                     !< Total number of dead pools

    !< Available pools
    REAL(wp), ALLOCATABLE :: dead(:,:)        !< (#nc, #dead pools): dead pools = woody plus livegrass for 1h deadfuel as proposed in Thonicke et al. (2010) [gDM/m2]
    REAL(wp), ALLOCATABLE :: live(:,:)        !< (#nc, #live pools): live pools= leaf, fruit, sap wood, heart wood(tree only), labile, reserve (grass and tree) [gDM/m2]
    REAL(wp), ALLOCATABLE :: live_frac(:,:)   !< Precompute fuel fraction allocation from each veg_pool_mt (#fuel classes, #live pools) [0 to 1], sum(live_frac,pool) = 1
    REAL(wp), ALLOCATABLE :: dead_frac(:,:)   !< Precompute fuel fraction allocation from each sb_pool_mt (#fuel classes, #dead pools) [0 to 1], sum(dead_frac,pool) = 1
  CONTAINS
    PROCEDURE :: Calc_dead_fuel               !< (#nc, #dead pools, #fuel classes) [g DM m-2]
!    PROCEDURE :: Calc_live_fuel               !< (#nc, #live pools, #fuel classes) [g DM m-2]
    PROCEDURE :: Calc_live_grass_fuel         !< (#nc, #live pools, #fuel classes) [g DM m-2]
    PROCEDURE :: Calc_live_tree_fuel          !< (#nc, #live pools, #fuel classes) [g DM m-2]

    !< Same as above but aggregated
!    PROCEDURE :: Calc_total_live_fuel         !< (#nc) [gDM m-2]
    PROCEDURE :: Calc_total_live_tree_fuel    !< (#nc) [gDM m-2]
    PROCEDURE :: Calc_total_live_grass_fuel   !< (#nc) [gDM m-2]

    PROCEDURE :: Count_dead_fuel_pools
  END TYPE t_sf_fuel_class


  INTERFACE t_sf_fuel_class
    MODULE PROCEDURE :: create_fuel_class
  END INTERFACE t_sf_fuel_class


CONTAINS

  ! ======================================================================================================= !
  !>
  !> Constructor for t_sf_fuel_class
  !>
  !> Either tree or grass is allocated according to "growth_form"
  !>
  FUNCTION create_fuel_class(nc, growth_form, dead_fuel_frac, live_fuel_frac, &
                                  dead_fuel_mt, live_fuel_mt, &
                                  litter_ratios, &
                                  cured_lvg_fract) &
                                  RESULT (fuel)
    INTEGER, INTENT(in)           :: nc
    INTEGER, INTENT(in)           :: growth_form(:)                  !< itree, igrass (#nc)
    REAL(wp), INTENT(in)          :: dead_fuel_frac(:,:)             !< Allocation of fuel from sb_pool_mt [0 to 1]
    REAL(wp), INTENT(in)          :: live_fuel_frac(:,:)             !< Allocation of fuel from veg_pool_mt [0 to 1]
    REAL(wp), INTENT(in)          :: live_fuel_mt(:,:)               !< bgcm live fuel (#pool, #nc) [g DM m-2]
    REAL(wp), INTENT(in)          :: dead_fuel_mt(:,:)               !< bgcm dead fuel (#pool, #nc) [g DM m-2]
    TYPE(t_sf_litter_ratios), INTENT(in) :: litter_ratios            !< Litter ratios for sb_pool [0 to 1]
    REAL(wp), INTENT(in)          :: cured_lvg_fract(:)              !< Cured live grass ratio [0 to 1] (#nc)
    TYPE(t_sf_fuel_class)         :: fuel                            !< output [gDM m-2]

    fuel%nc = nc
    ALLOCATE(fuel%growth_form, source = growth_form)
    fuel%litter_ratios = litter_ratios

    fuel%ndeadpools = sf_fuel_dead_fuel_nr !< woody litter + live grass pools (+ polymeric litter, won't be used if woody litter only)

    ALLOCATE(fuel%cured_lvg_fract(nc))
    fuel%cured_lvg_fract = cured_lvg_fract

    ALLOCATE(fuel%live(nc, sf_fuel_live_fuel_nr))
    ALLOCATE(fuel%dead(nc, fuel%ndeadpools))
    fuel%live = 0._wp
    fuel%dead = 0._wp

    ALLOCATE(fuel%dead_frac, source=dead_fuel_frac)
    ALLOCATE(fuel%live_frac, source=live_fuel_frac)

    !< Select sb_pool_mt pools as dead fuel
    fuel%dead(:,t_sf_fuel_dead_fuel_woody_ix) = dead_fuel_mt(ix_woody_litter, :)
    fuel%dead(:,t_sf_fuel_dead_fuel_polymeric_ix) = dead_fuel_mt(ix_polymeric_litter, :) * fuel%litter_ratios%polymeric
    fuel%dead(:,t_sf_fuel_dead_fuel_dom_assoc_ix) = dead_fuel_mt(ix_dom_assoc,:) * fuel%litter_ratios%dom_assoc
    fuel%dead(:,t_sf_fuel_dead_fuel_dom_ix) = dead_fuel_mt(ix_dom,:) * fuel%litter_ratios%dom
    fuel%dead(:,t_sf_fuel_dead_fuel_soluable_ix) = dead_fuel_mt(ix_soluable_litter,:) * fuel%litter_ratios%soluable
    fuel%dead(:,t_sf_fuel_dead_fuel_res_ix) = dead_fuel_mt(ix_residue,:) * fuel%litter_ratios%residue
    fuel%dead(:,t_sf_fuel_dead_fuel_res_assoc_ix) = dead_fuel_mt(ix_residue_assoc,:) * fuel%litter_ratios%residue_assoc

    !< Select veg_pool_mt pools as live fuel
    fuel%live(:,t_sf_fuel_live_fuel_leaf_ix) = live_fuel_mt(ix_leaf,:)
    fuel%live(:,t_sf_fuel_live_fuel_sap_wood_ix) = live_fuel_mt(ix_sap_wood,:)
    fuel%live(:,t_sf_fuel_live_fuel_heart_wood_ix) = live_fuel_mt(ix_heart_wood,:)
    fuel%live(:,t_sf_fuel_live_fuel_labile_ix) = live_fuel_mt(ix_labile,:)
    fuel%live(:,t_sf_fuel_live_fuel_reserve_ix) = live_fuel_mt(ix_reserve,:)
    fuel%live(:,t_sf_fuel_live_fuel_fruit_ix) = live_fuel_mt(ix_fruit,:)
  END FUNCTION create_fuel_class


  ! ======================================================================================================= !
  !>
  !> Calculate dead fuel by fuel class
  !> Following Thonicke et al. (2010)
  !> -Woody litter is divided into 1, 10, 100 and 1000 hours fuel
  !> -live grass is included into dead fuel only, all to 1 hours fuel
  !> -Polymeric litter is customized and experimental, so optional.
  !>
  FUNCTION Calc_dead_fuel(this) RESULT (dead_fuel)
    CLASS(t_sf_fuel_class), INTENT(in) :: this
    REAL(wp), ALLOCATABLE :: dead_fuel(:,:,:) !< (this%nc, this%ndeadpools, df_sf_fc_full_nr) !< output [g DM m-2]

    INTEGER :: ifuel

    IF (.NOT. ALLOCATED(dead_fuel)) ALLOCATE(dead_fuel(this%nc, this%ndeadpools, df_sf_fc_full_nr))

    dead_fuel = 0._wp

    DO ifuel=1, df_sf_fc_full_nr !< iterate over fuel classes (1,10,100,1000 hours)
      !< woody litter, polymeric, dom_assoc, dom, soluable litter, residue, residue_assoc
      dead_fuel(:,t_sf_fuel_dead_fuel_woody_ix,ifuel) = this%dead(:,t_sf_fuel_dead_fuel_woody_ix) * &
                              this%dead_frac(ifuel,t_sf_fuel_dead_fuel_woody_ix)
      dead_fuel(:,t_sf_fuel_dead_fuel_polymeric_ix,ifuel) = this%dead(:,t_sf_fuel_dead_fuel_polymeric_ix) * &
                              this%dead_frac(ifuel,t_sf_fuel_dead_fuel_polymeric_ix)
      dead_fuel(:,t_sf_fuel_dead_fuel_dom_assoc_ix,ifuel) = this%dead(:,t_sf_fuel_dead_fuel_dom_assoc_ix) * &
                              this%dead_frac(ifuel,t_sf_fuel_dead_fuel_dom_assoc_ix)
      dead_fuel(:,t_sf_fuel_dead_fuel_dom_ix,ifuel) = this%dead(:,t_sf_fuel_dead_fuel_dom_ix) * &
                              this%dead_frac(ifuel,t_sf_fuel_dead_fuel_dom_ix)
      dead_fuel(:,t_sf_fuel_dead_fuel_soluable_ix,ifuel) = this%dead(:,t_sf_fuel_dead_fuel_soluable_ix) * &
                              this%dead_frac(ifuel,t_sf_fuel_dead_fuel_soluable_ix)
      dead_fuel(:,t_sf_fuel_dead_fuel_res_ix,ifuel) = this%dead(:,t_sf_fuel_dead_fuel_res_ix) * &
                              this%dead_frac(ifuel,t_sf_fuel_dead_fuel_res_ix)
      dead_fuel(:,t_sf_fuel_dead_fuel_res_assoc_ix,ifuel) = this%dead(:,t_sf_fuel_dead_fuel_res_assoc_ix) * &
                              this%dead_frac(ifuel,t_sf_fuel_dead_fuel_res_assoc_ix)

      !< livegrass: leaf, sap wood, labile, reserve, fruit
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_leaf_ix,ifuel) = this%live(:,t_sf_fuel_live_fuel_leaf_ix) * &
                                                        this%live_frac(ifuel,t_sf_fuel_live_fuel_leaf_ix) * &
                                                        this%cured_lvg_fract
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_sap_ix,ifuel) = this%live(:,t_sf_fuel_live_fuel_sap_wood_ix) * &
                                                        this%live_frac(ifuel,t_sf_fuel_live_fuel_sap_wood_ix) * &
                                                        this%cured_lvg_fract
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_labile_ix,ifuel) = this%live(:,t_sf_fuel_live_fuel_labile_ix) * &
                                                        this%live_frac(ifuel,t_sf_fuel_live_fuel_labile_ix) * &
                                                        this%cured_lvg_fract
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_reserve_ix,ifuel) = this%live(:,t_sf_fuel_live_fuel_reserve_ix) * &
                                                        this%live_frac(ifuel,t_sf_fuel_live_fuel_reserve_ix) * &
                                                        this%cured_lvg_fract
      dead_fuel(:,t_sf_fuel_dead_fuel_grass_fruit_ix,ifuel) = this%live(:,t_sf_fuel_live_fuel_fruit_ix) * &
                                                        this%live_frac(ifuel,t_sf_fuel_live_fuel_fruit_ix) * &
                                                        this%cured_lvg_fract
    END DO

  END FUNCTION Calc_dead_fuel


  ! ======================================================================================================= !
  !>
  !> Calculate live fuel by fuel pool and fuel class
  !>
!  FUNCTION Calc_live_fuel(this) RESULT (live_fuel)
!    CLASS(t_sf_fuel_class), INTENT(in) :: this
!    REAL(wp) :: live_fuel(this%nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr)

!    live_fuel = this % Calc_live_tree_fuel()
!  END FUNCTION Calc_live_fuel


  ! ======================================================================================================= !
  !>
  !> Calculate live tree fuel, by pool and fuel class
  !>
  FUNCTION Calc_live_tree_fuel(this) RESULT (tree_fuel)
    CLASS(t_sf_fuel_class), INTENT(in) :: this
    REAL(wp) :: tree_fuel(this%nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr)

    INTEGER :: ifuel, ic

    tree_fuel = 0._wp !< igrass = 0

    DO ic=1, this%nc
      IF (this%growth_form(ic) == itree) THEN
        DO ifuel = 1, df_sf_fc_full_nr
          tree_fuel(ic,t_sf_fuel_live_fuel_leaf_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_leaf_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_leaf_ix)
          tree_fuel(ic,t_sf_fuel_live_fuel_sap_wood_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_sap_wood_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_sap_wood_ix)
          tree_fuel(ic,t_sf_fuel_live_fuel_heart_wood_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_heart_wood_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_heart_wood_ix)
          tree_fuel(ic,t_sf_fuel_live_fuel_labile_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_labile_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_labile_ix)
          tree_fuel(ic,t_sf_fuel_live_fuel_reserve_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_reserve_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_reserve_ix)
          tree_fuel(ic,t_sf_fuel_live_fuel_fruit_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_fruit_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_fruit_ix)
        END DO

      END IF
    END DO
  END FUNCTION Calc_live_tree_fuel


  ! ======================================================================================================= !
  !>
  !> Calculate live grass fuel by pool and fuel class
  !>
  FUNCTION Calc_live_grass_fuel(this) RESULT (grass_fuel)
    CLASS(t_sf_fuel_class), INTENT(in) :: this
    REAL(wp) :: grass_fuel(this%nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr)

    INTEGER :: ifuel, ic

    grass_fuel = 0._wp !< for itree = 0

    DO ic = 1, this % nc
      IF (this % growth_form(ic) == igrass) THEN
        DO ifuel = 1, df_sf_fc_full_nr
          grass_fuel(ic,t_sf_fuel_live_fuel_leaf_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_leaf_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_leaf_ix) * (1._wp - this%cured_lvg_fract(ic))
          grass_fuel(ic,t_sf_fuel_live_fuel_sap_wood_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_sap_wood_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_sap_wood_ix) * (1._wp - this%cured_lvg_fract(ic))
          grass_fuel(ic,t_sf_fuel_live_fuel_labile_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_labile_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_labile_ix) * (1._wp - this%cured_lvg_fract(ic))
          grass_fuel(ic,t_sf_fuel_live_fuel_reserve_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_reserve_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_reserve_ix) * (1._wp - this%cured_lvg_fract(ic))
          grass_fuel(ic,t_sf_fuel_live_fuel_fruit_ix, ifuel) = this%live(ic,t_sf_fuel_live_fuel_fruit_ix) * &
              this%live_frac(ifuel,t_sf_fuel_live_fuel_fruit_ix) * (1._wp - this%cured_lvg_fract(ic))
        END DO
      END IF
    END DO
  END FUNCTION Calc_live_grass_fuel


  ! ======================================================================================================= !
  !>
  !> Calculate total live tree fuel
  !>
  !> In case growth_form is grass, then it returns 0
  !>
!  FUNCTION Calc_total_live_fuel(this) RESULT (total_live_fuel)
 !   CLASS(t_sf_fuel_class), INTENT(in) :: this
 !   REAL(wp) :: total_live_fuel(this%nc)          !< output
 !   REAL(wp), ALLOCATABLE :: live_fuel(:,:,:)       !< (#nc, #live fuel pools, #fuel class) [gDM m-2]!

!    ALLOCATE(live_fuel(this%nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr))
!    live_fuel = this%Calc_live_fuel()
!    total_live_fuel = sf_calc_aggr_fuels(live_fuel)
!  END FUNCTION Calc_total_live_fuel


  ! ======================================================================================================= !
  !>
  !> Calculate total live fuel
  !>
  !> In case growth_form is grass, then it returns 0
  !>
  FUNCTION Calc_total_live_tree_fuel(this) RESULT (total_live_fuel)
    CLASS(t_sf_fuel_class), INTENT(in) :: this
    REAL(wp) :: total_live_fuel(this%nc)          !< output
    REAL(wp), ALLOCATABLE :: live_fuel(:,:,:)       !< (#nc, #live fuel pools, #fuel class) [gDM m-2]

    ALLOCATE(live_fuel(this%nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr))
    live_fuel = this%Calc_live_tree_fuel()
    total_live_fuel = sf_calc_aggr_fuels(live_fuel)
  END FUNCTION Calc_total_live_tree_fuel


  ! ======================================================================================================= !
  !>
  !> Calculate live grass fuel
  !>
  !> In case growth_form is tree, then it returns 0
  !>
  FUNCTION Calc_total_live_grass_fuel(this) RESULT (total_live_fuel)
    CLASS(t_sf_fuel_class), INTENT(in) :: this
    REAL(wp) :: total_live_fuel(this%nc)      !< ouput [gDM m-2]
    REAL(wp), ALLOCATABLE :: live_fuel(:,:,:)   !< (#nc, #live fuel pools, #fuel class) [gDM m-2]

    ALLOCATE(live_fuel(this%nc, sf_fuel_live_fuel_nr, df_sf_fc_full_nr))
    live_fuel = this%Calc_live_grass_fuel()
    total_live_fuel = sf_calc_aggr_fuels(live_fuel)
  END FUNCTION Calc_total_live_grass_fuel


  ! ======================================================================================================= !
  !>
  !> True when total_dead_fuel has any fuel
  !> Implementation for dead_fuel with 3 dimensions(:,:,:)
  !>
  !> Note: In fortran a generic function must implemente each type and rank, See INTERFACE sf_calc_has_any_fuel
  !> todo: rename dead_fuel to fuel (a more generic naming).
  !>
  FUNCTION sf_calc_has_any_fuel_r3(dead_fuel) RESULT (has_dead_fuel)
    USE mo_jsb_math_constants,   ONLY: eps8

    REAL(wp), INTENT(in)  :: dead_fuel(:,:,:)      !< [g DM m-2]
    LOGICAL, ALLOCATABLE  :: has_dead_fuel(:)      !< True if gridcell has any fuel (> 0)
    REAL(wp), ALLOCATABLE :: sum_of_dead_fuels(:)
    INTEGER :: nc

    nc =SIZE(dead_fuel, dim=1) !< 1=nc index
    IF (.NOT. ALLOCATED(has_dead_fuel)) ALLOCATE(has_dead_fuel(nc))

    sum_of_dead_fuels = sf_calc_aggr_fuels(dead_fuel)
    has_dead_fuel = (sum_of_dead_fuels > eps8)
  END FUNCTION sf_calc_has_any_fuel_r3

  ! ======================================================================================================= !
  !>
  !> True when total_dead_fuel has any fuel
  !> Implementation for dead_fuel with 2 dimensions(:,:)
  !>
  !> Note: In fortran a generic function must implemente each type and rank, See INTERFACE sf_calc_has_any_fuel
  !>
  FUNCTION sf_calc_has_any_fuel_r2(dead_fuel) RESULT (has_dead_fuel)
    USE mo_jsb_math_constants,   ONLY: eps8

    REAL(wp), INTENT(in)  :: dead_fuel(:,:)        !< [g DM m-2]
    LOGICAL, ALLOCATABLE  :: has_dead_fuel(:)     !< True if gridcell has any fuel (> 0)
    INTEGER :: nc

    nc = SIZE(dead_fuel, dim=1) !< 1 = nc index
    IF (.NOT. ALLOCATED(has_dead_fuel)) ALLOCATE(has_dead_fuel(nc))

    has_dead_fuel = (SUM(dead_fuel, dim=2) > eps8)
  END FUNCTION sf_calc_has_any_fuel_r2


  ! ======================================================================================================= !
  !>
  !> Sum of all fuels into nc
  !> fuel(nc, fuel pools, fuel classes) -> fuel(nc)
  !>
  IMPURE FUNCTION sf_calc_aggr_fuels(fuel, ncdim) RESULT (total_fuel)
    REAL(wp), INTENT(in)          :: fuel(:,:,:)    !< Dead fuel (#nc, #fuel pools, #fuel classes)
    INTEGER, INTENT(in), OPTIONAL :: ncdim          !< NC dimension index
    REAL(wp), ALLOCATABLE         :: total_fuel(:)  !< output (#nc)

    INTEGER, PARAMETER :: nc_idx_def = 1 !< Default NC dimension index
    INTEGER :: inc, nc, nc_idx

    nc_idx = nc_idx_def
    IF (PRESENT(ncdim)) nc_idx = ncdim

    nc = SIZE(fuel, dim=nc_idx)
    IF (.NOT. ALLOCATED(total_fuel)) ALLOCATE(total_fuel(nc))
    total_fuel(:) = 0._wp

    !< sum all dimensions from "fuel"
    SELECT CASE (nc_idx)
    CASE (1)
      DO CONCURRENT (inc = 1:nc)
        total_fuel(inc) = sum(fuel(inc,:,:))
      END DO
    CASE (3)
      DO CONCURRENT (inc = 1:nc)
        total_fuel(inc) = sum(fuel(:,:,inc))
      END DO
    CASE DEFAULT
      WRITE(*,*) "mo_q_dfire_sf_fuel_class.f90:: sf_calc_aggr_fuels:: NC dim =", nc_idx
      CALL finish("sf_calc_aggr_fuels","NCDIM argument is not implemented")
    END SELECT

  END FUNCTION sf_calc_aggr_fuels


  ! ======================================================================================================= !
  !>
  !> Count the number of dead fuel pools
  !>
  FUNCTION Count_dead_fuel_pools(this) RESULT (nfuelpools)
    CLASS(t_sf_fuel_class), INTENT(in) :: this
    INTEGER :: nfuelpools
    nfuelpools = this % ndeadpools
  END FUNCTION Count_dead_fuel_pools

#endif
END MODULE mo_q_dfire_sf_fuel_class
