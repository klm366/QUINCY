!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!>
!> @brief module to determine the fuel emissions based on the Pfeiffer et al. (2013) paper
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_sf_emissions_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

  USE mo_dfire_sf_adaptor,    ONLY: sf_convert_units
  USE mo_dfire_constants,     ONLY: df_sf_species_co2_idx, df_sf_species_co_idx, df_sf_species_ch4_idx
  USE mo_dfire_constants,     ONLY: df_sf_species_nh3_idx, df_sf_species_n2o_idx, df_sf_species_n2_idx

  USE mo_lnd_bgcm_idx, ONLY: ixC, ixC13, ixC14, ixN, ixN15, ixP
  USE mo_lnd_bgcm_idx, ONLY: ix_woody_litter !< sb_pool
  USE mo_lnd_bgcm_idx, ONLY: ix_fruit, ix_leaf, ix_reserve, ix_labile, ix_sap_wood, ix_heart_wood !< veg_pool
  USE mo_dfire_constants, ONLY: t_sf_fuel_dead_fuel_woody_ix, & !< dead fuel (litter)
                                  t_sf_fuel_dead_fuel_polymeric_ix, &
                                  t_sf_fuel_dead_fuel_soluable_ix, &
                                  t_sf_fuel_dead_fuel_res_ix, &
                                  t_sf_fuel_dead_fuel_res_assoc_ix, &
                                  t_sf_fuel_dead_fuel_dom_ix, &
                                  t_sf_fuel_dead_fuel_dom_assoc_ix

  USE mo_dfire_constants, ONLY: t_sf_fuel_dead_fuel_grass_leaf_ix, & !< dead fuel (livegrass)
                                  t_sf_fuel_dead_fuel_grass_sap_ix, &
                                  t_sf_fuel_dead_fuel_grass_labile_ix, &
                                  t_sf_fuel_dead_fuel_grass_reserve_ix, &
                                  t_sf_fuel_dead_fuel_grass_fruit_ix

  USE mo_dfire_constants, ONLY: t_sf_fuel_live_fuel_leaf_ix, t_sf_fuel_live_fuel_sap_wood_ix, & ! live fuel
                                  t_sf_fuel_live_fuel_heart_wood_ix, t_sf_fuel_live_fuel_labile_ix, &
                                  t_sf_fuel_live_fuel_reserve_ix, t_sf_fuel_live_fuel_fruit_ix

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_emissions, sf_calc_emissions_elements

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_sf_emissions_class'

CONTAINS

  !>
  !> Constructor
  !>
  !> Based on Eqs. A83, A85 and A86 from Pfeiffer et al. (2013)
  !> Further details:
  !> https://gitlab.dkrz.de/quincy-community/quincy-model-developers/qs-development/-/issues/113#peiffers-eqs-a83-to-a86
  !>
  FUNCTION sf_calc_emissions(nc, npools, nspecs, emission_factors, bb, mask) RESULT (emissions)
    INTEGER, INTENT(in)            :: nc
    INTEGER, INTENT(in)            :: npools !< Total number of pools
    INTEGER, INTENT(in)            :: nspecs !< total number of species
    REAL(wp), INTENT(in)           :: emission_factors(:) !< emissions factors (#nspecs) [molX/kgDM/m2]
    REAL(wp), INTENT(in)           :: bb(:,:) !< Burned biomass (#nc,#npools) [gDM m-2]
    LOGICAL, INTENT(in)            :: mask(:)   !< True = compute gridcell(s) (#nc)
    REAL(wp), ALLOCATABLE          :: emissions(:,:,:) !< output (#nc,#npools,#nspecs) [mol m-2]

    INTEGER :: ispecs, ipool, ic

    IF (.NOT. ALLOCATED(emissions)) ALLOCATE(emissions(nc,npools,nspecs))
    emissions = 0._wp

    ! Eq. A86 from Pfeiffer et al. (2013)
    DO ispecs=1, nspecs
      DO ipool=1, npools
        DO ic=1, nc
          IF (mask(ic)) THEN
            emissions(ic,ipool,ispecs) = emission_factors(ispecs) *  sf_convert_units("g", "kg", bb(ic,ipool))
          END IF
        END DO
      END DO
    END DO
  END FUNCTION sf_calc_emissions


  ! ======================================================================================================= !
  !>
  !> Calculate residue from emissions vs fuel consumption
  !>
  !> Carbon and Nitrogen nutrients are directly calculated from emissions.
  !> P is not yet calculated. A scientific discussion will be required on this topic. Currently remains in soil.
  !> Isotopes are derived from Nutrients at the same ratio for emitted elements.
  !> Residue can be directly calculated "FC - emitted elements"
  !> Any non mass conservation is capped to 0 but included into mass overflow for diagnostics.
  !>
  SUBROUTINE sf_calc_emissions_elements(nc, npools, is_type_of_fuel, &
      fuel_consumption, emissions, & !< in
      !residue, mass_overflow, & !< out
      emissions_el, & !< out
      mask)  !< optional

    INTEGER, INTENT(in)   :: nc, npools
    INTEGER, INTENT(in)   :: is_type_of_fuel         !< 1=dead fuel, 2=live fuel
    REAL(wp), INTENT(in)  :: fuel_consumption(:,:,:) !< [molX/m2] (#nc,#npools, #elements)
    REAL(wp), INTENT(in)  :: emissions(:,:,:)        !< [molX/m2] (#nc,#npools,#nspecs)
    REAL(wp), INTENT(out) :: emissions_el(:,:,:)     !< Emissions by elements [molX/m2] (#nc,#npools,#elements)
    LOGICAL, INTENT(in), OPTIONAL :: mask(:)         !< True=compute gridcell(s)

    INTEGER               :: ipool, ic
    LOGICAL               :: mask_opt(nc)

    mask_opt = .TRUE. !< default
    IF (PRESENT(mask)) THEN
      mask_opt = mask
    END IF
    emissions_el = 0._wp

    !< Back calculate emitted elements into C, N and P = 0
    DO ipool = 1, npools
      DO ic=1, SIZE(emissions, DIM=1)
        IF (mask(ic)) THEN
          emissions_el(ic,ipool,ixC) = emissions_el(ic,ipool,ixC) + emissions(ic,ipool,df_sf_species_co2_idx)  !< CO2
          emissions_el(ic,ipool,ixC) = emissions_el(ic,ipool,ixC) + emissions(ic,ipool,df_sf_species_co_idx)  !< CO
          emissions_el(ic,ipool,ixC) = emissions_el(ic,ipool,ixC) + emissions(ic,ipool,df_sf_species_ch4_idx) !< CH4

          emissions_el(ic,ipool,ixN) = emissions_el(ic,ipool,ixN) + emissions(ic,ipool,df_sf_species_nh3_idx) !< NH3
          emissions_el(ic,ipool,ixN) = emissions_el(ic,ipool,ixN) + emissions(ic,ipool,df_sf_species_n2o_idx) !< N2O
          emissions_el(ic,ipool,ixN) = emissions_el(ic,ipool,ixN) + emissions(ic,ipool,df_sf_species_n2_idx) !< N2
        END IF
      END DO

      !< Note: no P emissions
      emissions_el(:,ipool,ixP) = 0._wp
    END DO

    IF (is_type_of_fuel == 1) THEN !< 1 = IS_DEAD_FUEL
      !< Populate emissions_el isotopes at the same ratio as its nutrient
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_woody_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_woody_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_grass_fruit_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_grass_fruit_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_grass_leaf_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_grass_leaf_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_grass_labile_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_grass_labile_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_grass_reserve_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_grass_reserve_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_grass_sap_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_grass_sap_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_polymeric_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_polymeric_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_soluable_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_soluable_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_res_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_res_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_res_assoc_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_res_assoc_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_dom_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_dom_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_dead_fuel_dom_assoc_ix,:), &
                                                  emissions_el(:,t_sf_fuel_dead_fuel_dom_assoc_ix,:), mask_opt)
    ELSE IF (is_type_of_fuel == 2) THEN !< 2 = IS_LIVE_FUEL
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_live_fuel_leaf_ix,:), &
                                                  emissions_el(:,t_sf_fuel_live_fuel_leaf_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_live_fuel_sap_wood_ix,:), &
                                                  emissions_el(:,t_sf_fuel_live_fuel_sap_wood_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_live_fuel_heart_wood_ix,:), &
                                                  emissions_el(:,t_sf_fuel_live_fuel_heart_wood_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_live_fuel_labile_ix,:), &
                                                  emissions_el(:,t_sf_fuel_live_fuel_labile_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_live_fuel_reserve_ix,:), &
                                                  emissions_el(:,t_sf_fuel_live_fuel_reserve_ix,:), mask_opt)
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,t_sf_fuel_live_fuel_fruit_ix,:), &
                                                  emissions_el(:,t_sf_fuel_live_fuel_fruit_ix,:), mask_opt)
    ELSE
      ! stop!
    END IF
  END SUBROUTINE sf_calc_emissions_elements


  ! ======================================================================================================= !
  !>
  !> Calculate amount for each isotope (C13,C14,N15) at the same ratio as its nutrient
  !>
  SUBROUTINE calc_isotopes_same_ratio_as_nutrients(var, amount, mask)
    REAL(wp), INTENT(in) :: var(:,:) !< (#nc,#elements)
    REAL(wp), INTENT(inout) :: amount(:,:) !< (#nc,#elements)
    LOGICAL, INTENT(in) :: mask(:) !< TRUE=compute gridcell(s) (#nc)
    INTEGER :: ic

    amount(:,ixC13) = 0._wp
    amount(:,ixC14) = 0._wp
    amount(:,ixN15) = 0._wp
    DO ic=1, SIZE(mask)
      IF (mask(ic)) THEN
        IF (var(ic,ixC) > 0._wp) amount(ic,ixC13)= var(ic,ixC13) * amount(ic,ixC)/var(ic,ixC)
        IF (var(ic,ixC) > 0._wp) amount(ic,ixC14)= var(ic,ixC14) * amount(ic,ixC)/var(ic,ixC)
        IF (var(ic,ixN) > 0._wp) amount(ic,ixN15)= var(ic,ixN15) * amount(ic,ixN)/var(ic,ixN)
      END IF
    END DO

  END SUBROUTINE calc_isotopes_same_ratio_as_nutrients

#endif
END MODULE mo_dfire_sf_emissions_class

