!> QUINCY Calculate fuel consumption
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
!> @brief module to determine the dead fuel consumption
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_sf_fuel_consump_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message, finish

  USE mo_dfire_constants,     ONLY: t_sf_fuel_dead_fuel_woody_ix, t_sf_fuel_dead_fuel_polymeric_ix, &
      t_sf_fuel_dead_fuel_grass_leaf_ix, &
      t_sf_fuel_dead_fuel_grass_sap_ix, t_sf_fuel_dead_fuel_grass_labile_ix, &
      t_sf_fuel_dead_fuel_grass_reserve_ix, t_sf_fuel_dead_fuel_grass_fruit_ix
  USE mo_dfire_constants,     ONLY: t_sf_fuel_live_fuel_leaf_ix, t_sf_fuel_live_fuel_sap_wood_ix, &
      t_sf_fuel_live_fuel_heart_wood_ix, t_sf_fuel_live_fuel_labile_ix, &
      t_sf_fuel_live_fuel_reserve_ix, t_sf_fuel_live_fuel_fruit_ix
  USE mo_lnd_bgcm_idx,        ONLY: ix_woody_litter, ix_sap_wood, ix_heart_wood, ix_labile, ix_reserve, &
      ix_fruit, ix_polymeric_litter, ix_leaf

  USE mo_lnd_bgcm_idx,        ONLY: ixC
  USE mo_dfire_sf_adaptor,    ONLY: sf_convert_units

  !< Total number of spitfire live fuel pools
  USE mo_dfire_constants,     ONLY: t_sf_fuel_live_fuel_nr

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_dead_fuel_consumption, sf_calc_live_fuel_consumption, &
      sf_calc_dead_fuel_consumption_populate_elements, &
      sf_calc_live_fuel_consumption_populate_elements

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_sf_fuel_consump_class'

CONTAINS

  !>
  !> Calculate dead fuel consumption
  !>
  !> Derived from equations B1, B4 and B5 from Thonicke et al. (2010)
  !>
  FUNCTION sf_calc_dead_fuel_consumption(area_burned_frac, &
      fract_consumption, dead_fuel, mask) RESULT (dead_fc)

    REAL(wp), INTENT(in)                 :: area_burned_frac(:)       !< Fraction of the area burned (#nc) [0 to 1]
    REAL(wp), INTENT(in)                 :: fract_consumption(:, :)   !< Dead fraction of consumption (#nc,#fuel classes) [0 to 1]
    REAL(wp), INTENT(in)                 :: dead_fuel(:,:,:)          !< (#nc, #dead fuel pools, #fuel classes) [g DM m-2]
    LOGICAL, INTENT(in)                  :: mask(:)                   !< True = gridcells to compute
    REAL(wp), ALLOCATABLE                :: dead_fc(:,:,:)            !< output (nc, ndeadpools, nfuelclass) [g DM m-2]

    INTEGER  :: iclass, idpool

    IF (.NOT. ALLOCATED(dead_fc)) ALLOCATE(dead_fc, mold=dead_fuel)
    dead_fc = 0._wp

    DO iclass=1, SIZE(dead_fuel, DIM=3) !< 3 = dead fuel idx
      DO idpool = 1, SIZE(dead_fuel, DIM=2) !< 2 = dead fuel pool indices
        WHERE (mask)
          !< Equations B1, B4 and B5 from Thonicke et al. (2010)
          dead_fc(:,idpool,iclass) = fract_consumption(:,iclass) * dead_fuel(:,idpool,iclass) * area_burned_frac(:)
        END WHERE
      ENDDO
    END DO
  END FUNCTION sf_calc_dead_fuel_consumption


  ! ======================================================================================================= !
  !>
  !> Compute live fuel consumption
  !>
  !> Eq. Appendix B2 from Thonicke et al. (2010)
  !>
  FUNCTION sf_calc_live_fuel_consumption(nc, nlivepools, nfuelclass, &
      prop_crown_scorch, fract_crown_scorch, &
      area_burned_frac, live_fuel, mask) RESULT (live_fc)
    INTEGER, INTENT(in) :: nc, nlivepools, nfuelclass
    REAL(wp), INTENT(in) :: prop_crown_scorch(nc)               !< Proportion of the crown scorch by the fire [0 to 1]
    REAL(wp), INTENT(in) :: fract_crown_scorch(nfuelclass)      !< Fraction of the crown scorch consumed as live fuel for each fuel class [0 to 1]
    REAL(wp), INTENT(in) :: area_burned_frac(nc)                !< Fraction of the area burned by a fire [0 to 1]
    REAL(wp), INTENT(in) :: live_fuel(nc,nlivepools,nfuelclass) !< Live fuel [g DM m-2]
    LOGICAL, INTENT(in)  :: mask(nc)                            !< True = compute gridcell(s)
    REAL(wp) :: live_fc(nc,t_sf_fuel_live_fuel_nr)              !< output [g DM m-2]

    INTEGER :: iclass, ilvpool

    live_fc = 0._wp

    DO ilvpool = 1, nlivepools
      DO iclass=1, nfuelclass
        WHERE (mask)
          live_fc(:,ilvpool) = live_fc(:,ilvpool) + (live_fuel(:,ilvpool,iclass) * prop_crown_scorch(:) * &
              area_burned_frac * fract_crown_scorch(iclass))
        ENDWHERE
      END DO
    END DO
  END FUNCTION sf_calc_live_fuel_consumption


  ! ======================================================================================================= !
  !>
  !> Populate remaining elements as dead fuel consumption is computed in DM.
  !>
  !> Apply same C ratio from dead fuel consumption to sb_pool to the remaining elements
  !>
  FUNCTION sf_calc_dead_fuel_consumption_populate_elements(nc, ndeadpools, &
      nEMelems, &
      veg_pool_mt, sb_pool_mt, &
      fc_dead, mask) RESULT (fc_dead_els)

    INTEGER, INTENT(in) :: nc, ndeadpools, nEMelems
    REAL(wp), INTENT(in) :: veg_pool_mt(:,:,:)        !< [molX/m2] (#pool, #elements, #nc)
    REAL(wp), INTENT(in) :: sb_pool_mt(:,:,:)         !< [molX/m3] (#pool, #elements, #nc)
    REAL(wp), INTENT(in) :: fc_dead(:,:)              !< Dead fuel consumption in QS units [molX/m2] (#nc, #ndeadpools)
    LOGICAL, INTENT(in)  :: mask(:)                   !< (#nc)
    REAL(wp), ALLOCATABLE :: fc_dead_els(:,:,:)       !< (nc,ndeadpools,df_sf_elements_nr) !< Output [molX/m2]

    INTEGER :: ic

    IF (.NOT. ALLOCATED(fc_dead_els)) THEN
      ALLOCATE(fc_dead_els(nc,ndeadpools,nEMelems))
    END IF
    fc_dead_els = 0._wp

    !< woody litter
    !write(*,*) "calc_dead_fc_pop_elems = fc_dead (woody) =", fc_dead(:,t_sf_fuel_dead_fuel_woody_ix)
    !write(*,*) "calc_dead_fc_pop_elems = sb_pool(woody)=", sb_pool_mt(ix_woody_litter,:,1)
    DO ic=1, SIZE(veg_pool_mt, DIM=3) !< 3=nc dimension index

      IF (.NOT. mask(ic)) CYCLE ! skip gridcell if no fire

      ! woody litter
      IF (sb_pool_mt(ix_woody_litter,ixC,ic) > 0._wp) THEN
        fc_dead_els(ic,t_sf_fuel_dead_fuel_woody_ix,:) = sb_pool_mt(ix_woody_litter,:,ic) * &
            fc_dead(ic,t_sf_fuel_dead_fuel_woody_ix) / sb_pool_mt(ix_woody_litter,ixC,ic)
      ENDIF

      ! polymeric litter
      IF (sb_pool_mt(ix_polymeric_litter,ixC,ic) > 0._wp) THEN
        fc_dead_els(ic,t_sf_fuel_dead_fuel_polymeric_ix,:) = sb_pool_mt(ix_polymeric_litter,:,ic) * &
            fc_dead(ic,t_sf_fuel_dead_fuel_polymeric_ix) / sb_pool_mt(ix_polymeric_litter,ixC,ic)
      END IF

      ! leaf
      IF (veg_pool_mt(ix_leaf,ixC,ic) > 0._wp) THEN
        fc_dead_els(ic,t_sf_fuel_dead_fuel_grass_leaf_ix,:) = veg_pool_mt(ix_leaf,:,ic) * &
            fc_dead(ic,t_sf_fuel_dead_fuel_grass_leaf_ix) / veg_pool_mt(ix_leaf,ixC,ic)
      ENDIF

      !< sap wood
      IF (veg_pool_mt(ix_sap_wood,ixC,ic) > 0._wp) THEN
        fc_dead_els(ic,t_sf_fuel_dead_fuel_grass_sap_ix,:) = veg_pool_mt(ix_sap_wood,:,ic) * &
            fc_dead(ic,t_sf_fuel_dead_fuel_grass_sap_ix) / veg_pool_mt(ix_sap_wood,ixC,ic)
      ENDIF

      !< labile
      IF (veg_pool_mt(ix_labile,ixC,ic) > 0._wp) THEN
        fc_dead_els(ic,t_sf_fuel_dead_fuel_grass_labile_ix,:) = veg_pool_mt(ix_labile,:,ic) * &
            fc_dead(ic,t_sf_fuel_dead_fuel_grass_labile_ix) / veg_pool_mt(ix_labile,ixC,ic)
      ENDIF

      !< reserve
      IF (veg_pool_mt(ix_reserve,ixC,ic) > 0._wp) THEN
        fc_dead_els(ic,t_sf_fuel_dead_fuel_grass_reserve_ix,:) = veg_pool_mt(ix_reserve,:,ic) * &
            fc_dead(ic,t_sf_fuel_dead_fuel_grass_reserve_ix) / veg_pool_mt(ix_reserve,ixC,ic)
      ENDIF

      !< fruit
      IF (veg_pool_mt(ix_fruit,ixC,ic) > 0._wp) THEN
        fc_dead_els(ic,t_sf_fuel_dead_fuel_grass_fruit_ix,:) = veg_pool_mt(ix_fruit,:,ic) * &
            fc_dead(ic,t_sf_fuel_dead_fuel_grass_fruit_ix) / veg_pool_mt(ix_fruit,ixC,ic)
      ENDIF
    END DO

  END FUNCTION sf_calc_dead_fuel_consumption_populate_elements


  ! ======================================================================================================= !
  !>
  !> Populate the remaining elements (N,P) and isotopes (C13,C14,N15) based on:
  !> -Carbon ratio from consumed live fuel and veg_pool for each pool
  !>
  FUNCTION sf_calc_live_fuel_consumption_populate_elements(nc, nlivepools, &
      nEMelems, &
      veg_pool_mt, &
      fc_live, mask) RESULT (fc_live_els)

    INTEGER, INTENT(in) :: nc, nlivepools, nEMelems
    REAL(wp), INTENT(in) :: veg_pool_mt(:,:,:)        !< [molX/m2] (#pool, #elements, #nc)
    REAL(wp), INTENT(in) :: fc_live(:,:)              !< Live fuel consumption in QS units [molX/m2] (#nc, #nlivepools)
    LOGICAL, INTENT(in)  :: mask(:)                   !< (#nc)
    REAL(wp), ALLOCATABLE :: fc_live_els(:,:,:)    !< (nc,nlivepools,df_sf_elements_nr) !< Output [molX/m3]
    INTEGER               :: ic

    IF (.NOT. ALLOCATED(fc_live_els)) THEN
      ALLOCATE(fc_live_els(nc,nlivepools,nEMelems))
    END IF
    fc_live_els = 0._wp

    DO ic=1, SIZE(veg_pool_mt, DIM=3) !< 3=nc dimension index

      IF (.NOT. mask(ic)) CYCLE ! skip gridcell if no fire

      !< leaf
      IF (veg_pool_mt(ix_leaf,ixC,ic) > 0._wp) THEN
        fc_live_els(ic,t_sf_fuel_live_fuel_leaf_ix,:) = veg_pool_mt(ix_leaf,:,ic) * &
            fc_live(ic,t_sf_fuel_live_fuel_leaf_ix) / veg_pool_mt(ix_leaf,ixC,ic)
      ENDIF
      ! sap wood
      IF (veg_pool_mt(ix_sap_wood,ixC,ic) > 0._wp) THEN
        fc_live_els(ic,t_sf_fuel_live_fuel_sap_wood_ix,:) = veg_pool_mt(ix_sap_wood,:,ic) * &
            fc_live(ic,t_sf_fuel_live_fuel_sap_wood_ix) / veg_pool_mt(ix_sap_wood,ixC,ic)
      ENDIF
      ! heart wood
      IF (veg_pool_mt(ix_heart_wood,ixC,ic) > 0._wp) THEN
        fc_live_els(ic,t_sf_fuel_live_fuel_heart_wood_ix,:) = veg_pool_mt(ix_heart_wood,:,ic) * &
            fc_live(ic,t_sf_fuel_live_fuel_heart_wood_ix) / veg_pool_mt(ix_heart_wood,ixC,ic)
      ENDIF
      ! labile
      IF (veg_pool_mt(ix_labile,ixC,ic) > 0._wp) THEN
        fc_live_els(ic,t_sf_fuel_live_fuel_labile_ix,:) = veg_pool_mt(ix_labile,:,ic) * &
            fc_live(ic,t_sf_fuel_live_fuel_labile_ix) / veg_pool_mt(ix_labile,ixC,ic)
      ENDIF
      ! reserve
      IF (veg_pool_mt(ix_reserve,ixC,ic) > 0._wp) THEN
        fc_live_els(ic,t_sf_fuel_live_fuel_reserve_ix,:) = veg_pool_mt(ix_reserve,:,ic) * &
            fc_live(ic,t_sf_fuel_live_fuel_reserve_ix) / veg_pool_mt(ix_reserve,ixC,ic)
      ENDIF
      ! fruit
      IF (veg_pool_mt(ix_fruit,ixC,ic) > 0._wp) THEN
        fc_live_els(ic,t_sf_fuel_live_fuel_fruit_ix,:) = veg_pool_mt(ix_fruit,:,ic) * &
            fc_live(ic,t_sf_fuel_live_fuel_fruit_ix) / veg_pool_mt(ix_fruit,ixC,ic)
      ENDIF

    END DO

  END FUNCTION sf_calc_live_fuel_consumption_populate_elements

#endif
END MODULE mo_dfire_sf_fuel_consump_class
