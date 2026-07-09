!>
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
MODULE mo_q_dfire_sf_emissions_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

  USE mo_q_dfire_sf_adaptor,  ONLY: sf_convert_units
  USE mo_q_dfire_sf_constants, ONLY: df_sf_species_co2_idx, df_sf_species_co_idx, df_sf_species_ch4_idx
  USE mo_q_dfire_sf_constants, ONLY: df_sf_species_nh3_idx, df_sf_species_n2o_idx, df_sf_species_n2_idx

  USE mo_lnd_bgcm_idx, ONLY: ixC, ixC13, ixC14, ixN, ixN15, ixP

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sf_calc_emissions ! main subroutine
  PUBLIC :: sf_calc_emissions_species, sf_calc_emitted_elements, sf_calc_emitted_isotopes, sf_patch_emissions

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_emissions_class'

CONTAINS

  !>
  !> Constructor
  !>
  !> Based on Eqs. A83, A85 and A86 from Pfeiffer et al. (2013)
  !> Further details:
  !> https://gitlab.dkrz.de/quincy-community/quincy-model-developers/qs-development/-/issues/113#peiffers-eqs-a83-to-a86
  !>
  FUNCTION sf_calc_emissions_species(emission_factors, bb, mask) RESULT (emissions)
    REAL(wp), INTENT(in)           :: emission_factors(:) !< emissions factors (#nspecs) [molX/kgDM/m2]
    REAL(wp), INTENT(in)           :: bb(:,:)             !< Burned biomass (#nc,#npools) [kgDM m-2]
    REAL(wp), ALLOCATABLE          :: emissions(:,:,:)    !< output (#nc,#npools,#nspecs) [mol m-2]
    LOGICAL, INTENT(in)            :: mask(:)             !< True = compute gridcell(s) (#nc)

    INTEGER :: ispecs, ipool, ic
    INTEGER :: nc, npools, nspecs


    nc = SIZE(bb, DIM=1) ! 1=nc index
    npools = SIZE(bb, DIM=2) ! 2=pools index
    nspecs = SIZE(emission_factors)

    IF (.NOT. ALLOCATED(emissions)) ALLOCATE(emissions(nc,npools,nspecs))
    emissions = 0._wp

    ! Eq. A86 from Pfeiffer et al. (2013)
    DO ispecs=1, nspecs
      DO ipool=1, npools
        DO ic=1, nc
          IF (mask(ic)) THEN
            emissions(ic,ipool,ispecs) = emission_factors(ispecs) *  bb(ic,ipool)
          END IF
        END DO
      END DO
    END DO
  END FUNCTION sf_calc_emissions_species


  ! ======================================================================================================= !
  !>
  !> Calculate elements/nutrients from species emitted
  !>
  SUBROUTINE sf_calc_emitted_elements(emissions_spc, emissions_els, mask)
    REAL(wp), INTENT(in) :: emissions_spc(:,:,:)    !< [molX/m2] (#nc,#npools,#nspecs)
    REAL(wp), INTENT(inout) :: emissions_els(:,:,:) !< [molX/m2] (#nc,#npools,#nelements)
    LOGICAL, INTENT(in), OPTIONAL :: mask(:)

    LOGICAL, ALLOCATABLE :: mask_opt(:)
    INTEGER :: nc, npools
    INTEGER :: ic, ipool ! counter

    nc=SIZE(emissions_spc, DIM=1) ! 1=nc index
    npools=SIZE(emissions_spc, DIM=2) ! 2= npools index

    ALLOCATE(mask_opt(SIZE(emissions_spc, DIM=1))) ! 1=nc index
    mask_opt = .TRUE. !< default
    IF (PRESENT(mask)) THEN
      mask_opt = mask
    END IF

    emissions_els = 0._wp

    !< Back calculate emitted elements into C, N and P = 0
    DO ipool = 1, npools
      DO ic=1, nc
        IF (mask(ic)) THEN
          emissions_els(ic,ipool,ixC) = emissions_els(ic,ipool,ixC) + emissions_spc(ic,ipool,df_sf_species_co2_idx)  !< CO2
          emissions_els(ic,ipool,ixC) = emissions_els(ic,ipool,ixC) + emissions_spc(ic,ipool,df_sf_species_co_idx)  !< CO
          emissions_els(ic,ipool,ixC) = emissions_els(ic,ipool,ixC) + emissions_spc(ic,ipool,df_sf_species_ch4_idx) !< CH4

          emissions_els(ic,ipool,ixN) = emissions_els(ic,ipool,ixN) + emissions_spc(ic,ipool,df_sf_species_nh3_idx) !< NH3
          emissions_els(ic,ipool,ixN) = emissions_els(ic,ipool,ixN) + emissions_spc(ic,ipool,df_sf_species_n2o_idx) !< N2O
          emissions_els(ic,ipool,ixN) = emissions_els(ic,ipool,ixN) + emissions_spc(ic,ipool,df_sf_species_n2_idx) !< N2
        END IF
      END DO
    END DO
  END SUBROUTINE sf_calc_emitted_elements


  ! ======================================================================================================= !
  !>
  !> Calculate emitted isotopes from emissions
  !> Isotopes are derived at the same ratio as fuel consumption
  !>
  SUBROUTINE sf_calc_emitted_isotopes(fuel_consumption, emissions_els, mask)
    REAL(wp), INTENT(in)    :: fuel_consumption(:,:,:)  !< [molX/m2] (#nc,#npools, #elements)
    REAL(wp), INTENT(inout) :: emissions_els(:,:,:)     !< Emissions by elements [molX/m2] (#nc,#npools,#elements)
    LOGICAL, INTENT(in), OPTIONAL :: mask(:)            !< True=compute gridcell(s)

    INTEGER :: nc, npools
    INTEGER :: ipool
    LOGICAL, ALLOCATABLE :: mask_opt(:)

    nc=SIZE(fuel_consumption, DIM=1) ! 1=nc index
    npools=SIZE(fuel_consumption, DIM=2) ! 2=npools index

    ALLOCATE(mask_opt(nc))
    mask_opt=.TRUE. ! default
    IF (PRESENT(mask)) mask_opt=mask

    DO ipool=1, npools
      CALL calc_isotopes_same_ratio_as_nutrients(fuel_consumption(:,ipool,:), emissions_els(:,ipool,:), mask_opt)
    END DO
  END SUBROUTINE sf_calc_emitted_isotopes


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


  ! ======================================================================================================= !
  !>
  !> This is temporary patch on Emissions until further research is done on this matter.
  !> Currently, there is an excess of emissions due to the given Emissions factors.
  !> This results as emitted N is > than its fuel consumption.
  !> For instance, mass balance is not preserved.
  !>
  !> Current patch caps emissions to its max potential value, defined in Fuel consumption.
  !> However, it breaks CN ratio, as well as other potential issues.
  !>
  SUBROUTINE sf_patch_emissions(fc, emissions, mask)
    REAL(wp), INTENT(in) :: fc(:,:,:)           ! Fuel consumption (nc,pools,elements)
    REAL(wp), INTENT(inout) :: emissions(:,:,:) ! Emissions (nc,pools,elements)
    LOGICAL, INTENT(in) :: mask(:)              !

    INTEGER :: ic, ip, iel

    DO iel=1, SIZE(emissions, DIM=3)
      DO ip=1, SIZE(emissions, DIM=2)
        DO ic=1, SIZE( emissions, DIM=1)
          IF ( mask(ic) .AND. fc(ic,ip,iel) < emissions(ic,ip,iel)) THEN !< Mass imabalance found
            emissions(ic,ip,iel) = fc(ic,ip,iel)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE sf_patch_emissions


  ! ======================================================================================================= !
  !>
  !> Calc emissions (by species and elements)
  !> Carbon and Nitrogen nutrients are directly calculated from emissions.
  !> P is not yet calculated. A scientific discussion will be required on this topic. Currently remains in soil.
  !> Isotopes are derived from Nutrients at the same ratio for emitted elements.
  !>
  !> Any non mass conservation is capped to 0 but included into mass overflow for diagnostics.
  !>
  SUBROUTINE sf_calc_emissions(emission_factors, bb, fc_els, emissions_spc, emissions_els, mask)
    REAL(wp), INTENT(in)          :: emission_factors(:)      !< emissions factors (#nspecs) [molX/kgDM/m2]
    REAL(wp), INTENT(in)          :: bb(:,:)                  !< Burnt biomass (#nc,#npools) [gDM m-2]
    REAL(wp), INTENT(in)          :: fc_els(:,:,:)            !< Fuel consumption [molX/m2] (#nc,#npools, #elements)
    REAL(wp), INTENT(out)         :: emissions_spc(:,:,:)     !< Emissions by species [molX/m2] (#nc,#npools,#elements)
    REAL(wp), INTENT(out)         :: emissions_els(:,:,:)     !< Emissions by elements [molX/m2] (#nc,#npools,#elements)
    LOGICAL, INTENT(in), OPTIONAL :: mask(:)                  !< True=compute gridcell(s)

    REAL(wp), ALLOCATABLE         :: bb_in_kg(:,:)            !< Burnt biomass (#nc,#npools) [kgDM m-2]

    ! Calc Emissions (by Species)
    !< Calculate fire damage to plants (total biomass burned and tree mortality) based on total dead and
    !< live fuel, as in Peiffer et al. (2013)
    bb_in_kg = sf_convert_units("g", "kg", bb)
    emissions_spc = sf_calc_emissions_species(emission_factors, bb_in_kg, mask)

    !< Convert emissions from species to elements (CNP), no emissions for P
    CALL sf_calc_emitted_elements(emissions_spc, emissions_els, mask)

    !< Account for isotopes, the same ratio found in fuel consumption
    CALL sf_calc_emitted_isotopes(fc_els, emissions_els, mask)

    !< Cap emissions to FC due to excess of matter in emissions.
    !< This avoids negative pools until it is solved.
    CALL sf_patch_emissions(fc_els, emissions_els, mask)
  END SUBROUTINE sf_calc_emissions

#endif
END MODULE mo_q_dfire_sf_emissions_class
