!> Contains the routines for the carbon processes
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
MODULE mo_carbon_process
#ifndef __NO_JSBACH__

  USE mo_jsb_control,      ONLY: acc_stream
  USE mo_kind,             ONLY: lp
  USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround, &
    & LCT_ACID, LCT_WATER, LCT_ETHANOL,                &
    & LCT_NONSOLUBLE, LCT_HUMUS,                       &
    & molarMassCO2_kg

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_plant_c, calc_soil_c, relocate_carbon_fire, relocate_carbon_damage, calc_vertical_soil_c, &
    & relocate_carbon_vertical_fire, relocate_carbon_vertical_damage, determine_root_dist, &
    & get_per_tile, add_litter_to_yasso_pool

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_process'

CONTAINS

  ! Calculations to allocate NPP to the plant C pools as well as various respiration and litter pools.
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_plant_c( &
    ! Inputs
    & delta_time_C, is_crop,      &
    & LAI,                        &
    & LAI_previous,               &
    & max_LAI,                    &
    & NPP_pot_yDayMean,           &
    !
    & fract_npp_2_woodPool,       &
    & fract_npp_2_reservePool,    &
    & fract_npp_2_exudates,       &
    & fract_green_2_herbivory,    &
    & tau_c_woods,                &
    & LAI_shed_constant,          &
    & Max_C_content_woods,        &
    & specific_leaf_area_C,       &
    & reserveC2leafC,             &
    ! Inouts
    & c_green,                    &
    & c_woods,                    &
    & c_reserve,                  &
    & c_crop_harvest,             &
    ! Outputs
    & NPP_flux_correction,        &
    & NPP_excess,                 &
    & NPP_act,                    &
    & Cflx_2_root_exudates,       &
    & Cflx_greenwood_2_litter,    &
    & Cflx_green_2_litter,        &
    & Cflx_wood_2_litter,         &
    & Cflx_herbivory,             &
    & Cflx_herb_2_littergreen,    &
    & Cflx_herb_2_atm,            &
    & Cflx_crop_harvest_2_atm     &
    & )
    !$ACC ROUTINE SEQ

    USE mo_carbon_constants, ONLY: fract_green_aboveGround
    !

    !! Names of carbon variables and arguments are:
    !!   1) Cflx_*  : Carbon flux [mol(C) m-2 s-1]
    !!   2) NPP*    : as 1)
    !!   3) c_*_2_* : Carbon transport [mol(C) m-2 carbon-time-step-1], where carbon-time-step is normally 1 day
    !!   4) c_*     : (if not covered by 3)) Carbon pool [mol(C) m-2]
    !! All [m-2] are refere to canopy area

    ! Input arguments
    REAL(lp), INTENT(in)    :: delta_time_C            !< Carbon time step [s]
    LOGICAL,  INTENT(in)    :: is_crop                 !< Is this a crop tile?
    REAL(lp), INTENT(in)    :: LAI                     !< Yesterdays mean LAI
    REAL(lp), INTENT(in)    :: LAI_previous            !< The day before yesterdays mean LAI
    REAL(lp), INTENT(in)    :: max_LAI                 !< Maximum value of LAI
    REAL(lp), INTENT(in)    :: NPP_pot_yDayMean        !< Yesterdays mean NPP rate
    REAL(lp), INTENT(in)    :: fract_npp_2_woodPool    !< Fraction of NPP to be put maximally into the green pool
    REAL(lp), INTENT(in)    :: fract_npp_2_reservePool !< Fraction of NPP to be put into the optimally into the reserve pool
    REAL(lp), INTENT(in)    :: fract_npp_2_exudates    !< Fraction of NPP to be put into the optimally into the root exudates
    REAL(lp), INTENT(in)    :: fract_green_2_herbivory
    REAL(lp), INTENT(in)    :: tau_c_woods             !< Time constant by which woods Pool is depreciated  [d]
                                                       !< Note: in lctlib it is given in [yr-1] but used by JS4 in [d-1].
    REAL(lp), INTENT(in)    :: LAI_shed_constant       !< Leaf shedding at a constant rate for evergreens etc. [d-1]
    REAL(lp), INTENT(in)    :: Max_C_content_woods     !< Maximum carbon content of wood pool (from lctLib)
    REAL(lp), INTENT(in)    :: specific_leaf_area_C    !< Specific leaf area (from lctLib) [m2(leaf) mol(C)-1]
    REAL(lp), INTENT(in)    :: reserveC2leafC          !< Ratio of max. carbon in reserve pool to max. carbon of leaves
    ! Plant carbon pools: on input: last value; updated on output
    REAL(lp), INTENT(inout) :: c_green                 !< Green carbon pool
    REAL(lp), INTENT(inout) :: c_woods                 !< Wood carbon pool
    REAL(lp), INTENT(inout) :: c_reserve               !< Reserve carbon pool
    REAL(lp), INTENT(inout) :: c_crop_harvest          !< Crop harvest carbon pool
    ! Output arguments
    REAL(lp), INTENT(out)   :: NPP_flux_correction     !< Amount by which the NPP rate entering the routine has to be corrected. This
                                                       !<   correction arises either because otherwise the reserve pool would get
                                                       !<   negative (positive correction), or the wood pool would exceed its
                                                       !<   maximum value (negative correction).
    REAL(lp), INTENT(out)   :: NPP_excess              !< Part of NPP that could not be stored in one of the plant carbon pools
                                                       !<   (green, wood, reserve), but had to be  thrown away (a posteriori
                                                       !<   reduction of NPP)
    REAL(lp), INTENT(out)   :: NPP_act                 !< Actual NPP after N-limitation and excess carbon drop
    REAL(lp), INTENT(out)   :: Cflx_2_root_exudates    !< Carbon flux into root exudates
    REAL(lp), INTENT(out)   :: Cflx_greenwood_2_litter !< Total carbon flux from vegetation to litter (green+woody, ag+bg)
    REAL(lp), INTENT(out)   :: Cflx_green_2_litter     !< Total carbon flux of easily decomposable material
    REAL(lp), INTENT(out)   :: Cflx_wood_2_litter      !< Total carbon flux of slowly decomposable material
    REAL(lp), INTENT(out)   :: Cflx_herbivory
    REAL(lp), INTENT(out)   :: Cflx_herb_2_littergreen
    REAL(lp), INTENT(out)   :: Cflx_herb_2_atm
    REAL(lp), INTENT(out)   :: Cflx_crop_harvest_2_atm

    ! Parameters
    REAL(lp), PARAMETER :: days_per_year      = 365.0_lp   !< Just for testing - should be imported, but slightly different value
    ! R: Normally parameters should be in mo_carbon_constants, but the following parameters are used only in this subroutine!
    !    Therefore I put them here. Is there any not structural advantage to put them into mo_carbon_constants?
    REAL(lp), PARAMETER :: greenC2leafC         = 4.0_lp          !< ratio of carbon green pool (leaves, fine roots, starches,
                                                                  !<   sugars) to leaf carbon
    REAL(lp), PARAMETER :: fract_c_crop_harvest = 0.5_lp          !< fraction of C from ag litter flux put to the crop harvest pool
    REAL(lp), PARAMETER :: fract_c_faeces2_LG   = 0.3_lp          !< fraction of C from faces put into the litter green pool
    REAL(lp), PARAMETER :: tau_c_reserve        = 1.0_lp * days_per_year !< time constant by which reserve pool is depreciated [yr]
    REAL(lp), PARAMETER :: tau_c_crop_harvest   = 1.0_lp * days_per_year !< time constant by which the crop harvest pool decays [yr]

    REAL(lp) :: c_green_max, c_reserve_optimal
    REAL(lp) :: c_NPP_2_greenPool, c_NPP_2_woodPool, c_NPP_2_reservePool, c_NPP_2_rootExudates, c_green_2_herbivory
    REAL(lp) :: c_2_litter_greenPool, c_2_litter_woodPool, c_reserve_2_litter_greenPool
    REAL(lp) :: c_excess, c_crop_harvest_2_atm
    REAL(lp) :: leaf_shedding

    REAL(lp) :: c_green_pot, c_woods_pot         !! values of potential carbon pools i.e. before accounting for N-limitation
    REAL(lp) :: c_reserve_pot

    REAL(lp) :: c_NPP_2_green_pot, c_NPP_2_wood_pot, c_NPP_2_reserve_pot
    REAL(lp) :: c_2_crop_harvest
    REAL(lp) :: c_faeces_2_LG, c_faeces_2_atm
    ! Preparation for N-limitation
    REAL(lp) :: c_NPP_2_green_via_mobileN, c_NPP_2_wood_via_mobileN
    REAL(lp) :: redFact_Nlimitation
    ! R: only for spitfire
    !REAL(lp) :: litter_wood_new, litter_wood_total

    ! Initializations
    c_2_crop_harvest     = 0.0_lp

    c_green_max       = greenC2leafC * LAI / specific_leaf_area_C
    c_reserve_optimal = reserveC2leafC * max_LAI / specific_leaf_area_C

    !! ============== NON-N-LIMITED CARBON AND NITROGEN ALLOCATION =========================================================
    !!
    !! Perform those C- and N-allocation steps that are not restricted by N-availability
    !! 1. Leaf shedding
    !! 2. Wood shedding
    !! 3. Depletion of the reserve pool
    !! 4. Decomposition of slow soil pool (whether this is really independent of N-availability is debatable!!)
    !!
    !!
    !! ----------------------------------------------------------------------------
    !! 1. Leaf shedding: Transfer of C from green pool to leaf litter litter pool

!!$ tr
!!$       IF (LAI >= LAI_previous) THEN                              !! If LAI is increasing or constant the leaf shedding is given
!!$          leaf_shedding = LAI * LAI_shed_constant                 !!    by a constant loss of leaves for evergreens, raingreens
!!$                                                                  !!    and grasses.
!!$       ELSE                                                       !! Otherwise
!!$          leaf_shedding = LAI_previous - LAI                      !!    the leaf shedding is given by the decrease in LAI.
!!$       END IF
!!$ tr

    leaf_shedding = MAX(LAI * LAI_shed_constant,LAI_previous - LAI) !! Leaf shedding is assured at a minimum constant rate
                                                                    !! which maybe exceeded by the decrease in LAI.
    c_2_litter_greenPool =  &                                       !! At most the whole carbon from the green pool is transferred
      & MIN(greenC2leafC * leaf_shedding / specific_leaf_area_c, &  !!    by leaf and root shedding to the green litter pool
      & c_green)
    c_green = c_green - c_2_litter_greenPool                        !! This removes carbon from the green pool
    c_excess = MAX(0.0_lp,c_green - c_green_max)                    !! If by the reduction of LAI the maximum value of the green
                                                                    !!    pool is still smaller than the current value, there is
    c_2_litter_greenPool = c_2_litter_greenPool + c_excess          !!    excess C that has also to be shedded to the green litter
    c_green              = c_green              - c_excess          !!    pool and also subtracted from the green pool

    IF (is_crop) THEN
      !!  Move shedded Carbon pro rata to crop harvest flux
      c_2_crop_harvest = fract_green_aboveGround * fract_c_crop_harvest * c_2_litter_greenPool
    END IF

    !! ---------------------------------------------------------------------------------
    !! Herbivory loss from c_green & Npool_green w.r.t PFTs in lctlib

    c_green_2_herbivory = c_green * fract_green_2_herbivory                   !! Grazing flux indirect comput.from NPP losses
    c_green             = c_green - c_green_2_herbivory                       !! Loss due to grazing

    c_faeces_2_LG       = c_green_2_herbivory * fract_c_faeces2_LG            !! c_faeces_2_LG to be stored in litter green pool
                                                                              !!   (not considered in Cflx_greenwood_2_litter)
                                                                              !!   (no life times (tau) enter c_faeces_2_LG!)
    c_faeces_2_atm      = c_green_2_herbivory * (1.0_lp - fract_c_faeces2_LG) !! c_faeces_2_atm to be lost to atmosphere
                                                                              !!   (dk: not added to CO2?)
    !! ---------------------------------------------------------------------------------
    !! 2. Wood shedding: Transfer of C from wood pool to wood litter pools
    !!
                                                            !! Assuming that forests continuously die at the inverse lifetime of
                                                            !! trees
    c_2_litter_woodPool = c_woods / tau_c_woods             ! The bigger tau gets, the weaker the C pools will be depreciated:
                                                            ! tau_c_woods: [days] -> [timestep]
                                                            !! .. the shedded wood (MAX FUNCTION FOR NONWOODY PFTS)
    ! R: only for spitfire
    ! if (present(fract_litter_wood_new)) litter_wood_new = fract_wood_aboveGround * c_2_litter_woodPool
    c_woods = c_woods -  c_2_litter_woodPool                !! .. and then all is subtracted from the wood pool.

    !! ----------------------------------------------------------------------------
    !! 3. Depletion of reserve pool: Transfer of C from reserve pool to green litter pool
    !!    Some organic carbon of the reserve pool is always lost with mortality of plants.
    !!    Note that the reserve pool contains no Nitrogen (starches and sugar are free of N).
    !!    Note also that the green litter pool has no fixed C/N-ratio (therefore no compensation
    !!    flux is needed)

    c_reserve_2_litter_greenPool = c_reserve / tau_c_reserve
    IF (is_crop) THEN
      c_2_crop_harvest = c_2_crop_harvest + fract_green_aboveGround * fract_C_crop_harvest * c_reserve_2_litter_greenPool
    END IF
    c_reserve      = c_reserve      - c_reserve_2_litter_greenPool
    c_crop_harvest = c_crop_harvest + c_2_crop_harvest

    !!
    !! ============== START OF POTENTIAL PLANT-CARBON ALLOCATION TO ESTIMATE PLANT-N-DEMAND =================================
    !!                       (exception: case NPP<0 handles actual allocation)

    c_green_pot   = c_green
    c_woods_pot   = c_woods
    c_reserve_pot = c_reserve

    !! Preliminary determination of distribution of NPP to the various pools (may be corrected afterwards according to
    !!    avaialable Nitrogen)

    NPP_flux_correction = 0.0_lp
    IF (NPP_pot_yDayMean <= 0.0_lp) THEN                      !! In case of negative or zero NPP
      NPP_excess = 0.0_lp

      c_NPP_2_reservePool =  NPP_pot_yDayMean * delta_time_C  !! .. the plant looses C and
      c_reserve     = c_reserve + c_NPP_2_reservePool         !! .. we try to take it from the reserve pool
      IF (c_reserve < 0.0_lp) THEN                            !! But if thereby the reserve pool gets negative
        c_NPP_2_reservePool = c_NPP_2_reservePool - c_reserve !! .. we absorb negative NPP only up to zero reserve pool and
        NPP_flux_correction = -c_reserve / delta_time_C       !! .. give the deficit as a flux correction back to the
        c_reserve = 0.0_lp                                    !! .. calling routine. Hence all Carbon from the reserve pool
      END IF                                                  !! .. is used up.

      !! For correct handling of nitrogen limitation below (which handles only the case of positive NPP) set NPP to zero:
      c_NPP_2_reservePool  = 0.0_lp
      c_NPP_2_greenPool    = 0.0_lp    !! All other transfer rates are zero
      c_NPP_2_woodPool     = 0.0_lp
      Cflx_2_root_exudates = 0.0_lp
      c_reserve_pot        = c_reserve !! remember modified reserve pool

    ELSE  !! NPP is positive

      !! NPP it is distributed according to predefined relative fraction
      c_NPP_2_woodPool     = delta_time_C * fract_npp_2_woodPool    * NPP_pot_yDayMean
      c_NPP_2_reservePool  = delta_time_C * fract_npp_2_reservePool * NPP_pot_yDayMean
      c_NPP_2_rootExudates = delta_time_C * fract_npp_2_exudates    * NPP_pot_yDayMean
      c_NPP_2_greenPool    = delta_time_C * &
        & (1.0_lp-fract_npp_2_woodPool-fract_npp_2_reservePool-fract_npp_2_exudates) * NPP_pot_yDayMean

      !! Growth of Wood Pool (structural carbon of living plants)

      c_woods_pot       = c_woods_pot + c_NPP_2_woodPool                !! Then it is attempted to put all NPP share into it
      c_excess          = MAX(0.0_lp,c_woods_pot-Max_C_content_woods)   !! .. but thereby the pool may get too large
      c_woods_pot       = c_woods_pot - c_excess                        !! .. so that the excess carbon has once more to be
                                                                        !! .. subtracted
      c_NPP_2_greenPool = c_NPP_2_greenPool + c_excess                  !! .. and instead made available to the green pool
      c_NPP_2_woodPool  = c_NPP_2_woodPool  - c_excess                  !! .. and the actual amount of NPP put to the wood pool
                                                                        !! .. is less

      !! Growth of Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
      IF (c_reserve_pot < c_reserve_optimal) THEN                       !! .. If the reserve pool is smaller than optimal,
        c_reserve_pot     = c_reserve_pot     + c_NPP_2_reservePool     !! .... it is filled by the available NPP
        c_excess          =                                         &   !! .... Thereby it may happen that it gets larger than
          & MAX(0.0_lp,c_reserve_pot - c_reserve_optimal)               !! .... optimal so that there is excess carbon
        c_reserve_pot       = c_reserve_pot       - c_excess            !! .... that needs not be taken up
        c_NPP_2_greenPool   = c_NPP_2_greenPool   + c_excess            !! .... but can better be used to increase the green
        c_NPP_2_reservePool = c_NPP_2_reservePool - c_excess            !! .... pool and the actual amount of NPP put to the
                                                                        !! .... reserve pool is less.
      ELSE                                                              !! .... Otherwise (reserve pool is larger than optimal)
        c_NPP_2_greenPool   = c_NPP_2_greenPool   + c_NPP_2_reservePool !! .... all NPP is left for the green pool
        c_NPP_2_reservePool = 0.0_lp                                    !! .... so that nothing is stored in the reserve pool.
      END IF

      !! Growth of Green Pool (leaves and fine roots): (in case of too much NPP, try to put it into the reserve pool).

      c_green_pot       = c_green_pot + c_NPP_2_greenPool               !! Green pool is filled by the available NPP.
      c_excess          = MAX(0.0_lp,c_green_pot - c_green_max)         !! .. Thereby it may get larger than appropriate for
                                                                        !!    current LAI.
      c_NPP_2_greenPool = c_NPP_2_greenPool - c_excess                  !! .. Hence the actual amount of NPP put to the green
                                                                        !!    pool is less.
      c_green_pot       = c_green_pot - c_excess                        !! .. This excess carbon needs not be taken up but
      IF (c_reserve_pot < c_reserve_optimal) THEN                       !! .... if the reserve pool is smaller than optimal
        c_reserve_pot       = c_reserve_pot       + c_excess            !! .... it is tried to put the carbon there,
        c_NPP_2_reservePool = c_NPP_2_reservePool + c_excess            !! .... which means that additional NPP is put to the
                                                                        !!      reserve pool.
        c_excess            =                                   &       !! .... Thereby it may happen that the reserve pool
          & MAX(0.0_lp,c_reserve_pot         - c_reserve_optimal)       !!      increases beyond the optimal value.
        c_reserve_pot       = c_reserve_pot       - c_excess            !! .... In that case the excess carbon is once more
                                                                        !!      removed from the reserve pool
        c_NPP_2_reservePool = c_NPP_2_reservePool - c_excess            !! .... so that the actual amount of NPP put into the
                                                                        !!      reserve pool is less.
      END IF

      NPP_excess           = c_excess             / delta_time_C
      Cflx_2_root_exudates = c_NPP_2_rootExudates / delta_time_C

    END IF !! NPP > 0 end

    !! ============== START OF POTENTIAL LITTER-CARBON ALLOCATION ==============================================

    !! ============== use nitrogen from N mobile pool to satisfy (at least partly) for potential carbon fluxes ===
    ! Remarks: Nmobile_2_green + Nmobile_2_wood <= Npool_mobile
    ! Npool_mobile accounts minN_plant_demand by trasfering N into green and wood pool

    c_NPP_2_green_pot = c_NPP_2_greenPool    !! potential growth of green carbon
    c_NPP_2_wood_pot  = c_NPP_2_woodPool     !! potential growth of wood  carbon

    c_NPP_2_green_via_mobileN = 0.0_lp
    c_NPP_2_wood_via_mobileN  = 0.0_lp

    c_NPP_2_green_pot = c_NPP_2_green_pot - c_NPP_2_green_via_mobileN
    c_NPP_2_wood_pot  = c_NPP_2_wood_pot  - c_NPP_2_wood_via_mobileN

    ! R: These lines were calculated also with_Nitrogen=.FALSE. in JS3. However, it was not used further when with_Nitrogen=.FALSE.
    !    In JS3 it should stand within the brackets.
    !Nflux_NPP_2_green_pot = c_NPP_2_green_pot/cn_green
    !Nflux_NPP_2_wood_pot  = c_NPP_2_wood_pot/cn_woods
    !minN_plant_demand = Nflux_NPP_2_green_pot + Nflux_NPP_2_wood_pot

    !! ============== COLLECT POTENTIAL C-FLUXES ===============================================================

    c_NPP_2_reserve_pot = c_NPP_2_reservePool !! potential growth of reserve carbon
    redFact_Nlimitation = 1.0_lp              !! This assures absence of nitrogen limitation in this mixed C-N code

    !! ============== UPDATE ALL C-POOLS =======================================================================

    c_green   = c_green   + c_NPP_2_green_pot * redFact_Nlimitation
    c_woods   = c_woods   + c_NPP_2_wood_pot  * redFact_Nlimitation
    c_reserve = c_reserve + c_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N

    !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

    IF (NPP_pot_yDayMean > 0.0_lp) THEN
      NPP_act = (c_NPP_2_green_via_mobileN + c_NPP_2_wood_via_mobileN + c_NPP_2_rootExudates + &
        & redFact_Nlimitation * (c_NPP_2_green_pot + c_NPP_2_wood_pot) + c_NPP_2_reserve_pot   &
        & ) / delta_time_C
    ELSE
      NPP_act = NPP_pot_yDayMean + NPP_flux_correction  ! NPP_pot_yDayMean is negative: NPP_flux_correction is zero
                                                         ! except c_reserve was < 0.0_lp.
    END IF

    !! ============== Decay of crop harvest =============

    c_crop_harvest_2_atm    =  c_crop_harvest / tau_c_crop_harvest
    c_crop_harvest          =  c_crop_harvest - c_crop_harvest_2_atm

    Cflx_crop_harvest_2_atm = -c_crop_harvest_2_atm / delta_time_C

    !! ============== COMPUTE Herbivory RESPIRATION =====

    Cflx_herb_2_atm = -(c_faeces_2_atm) / delta_time_C

    !! ============== Litter and herbivory diagnostics (flux from d-1 to s-1)
    Cflx_greenwood_2_litter = (c_2_litter_greenPool + c_2_litter_woodPool + c_reserve_2_litter_greenPool) / delta_time_C
    Cflx_green_2_litter     = Cflx_greenwood_2_litter + (-c_2_litter_woodPool+c_faeces_2_LG-c_2_crop_harvest) / delta_time_C
    Cflx_wood_2_litter      = c_2_litter_woodPool / delta_time_C
    Cflx_herbivory          = c_Green_2_herbivory / delta_time_C
    Cflx_herb_2_littergreen = c_faeces_2_LG  / delta_time_C   ! DSG: goes into yasso

    ! R: only for spitfire:
    ! fraction of new wood litter, assume old and new litter are respired in the same way, respiration does not change the fraction
    !        if (present(fract_litter_wood_new)) then
    ! litter_wood_total = c_water_ag2 + c_acid_ag2 + c_ethanol_ag2 + c_nonsoluble_ag2
    !           if ((litter_wood_total >  EPSILON(1._lp)) .AND. (litter_wood_new < litter_wood_total)) then
    !              fract_litter_wood_new = litter_wood_new / litter_wood_total
    !           else
    !              fract_litter_wood_new = 1._lp    ! fract_litter_wood_new=1 means fuel fractions will be reset to original values
    !           endif
    !        endif


  END SUBROUTINE calc_plant_c


  ! Calculations to allocate NPP to the plant C pools.
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE calc_soil_c( &
    ! General inputs
    & delta_time_C,               &
    ! Input coefficients (fixed)
    & LeafLit_coef_acid,          &
    & LeafLit_coef_water,         &
    & LeafLit_coef_ethanol,       &
    & LeafLit_coef_nonsoluble,    &
    & LeafLit_coef_humus,         &
    & WoodLit_coef_acid,          &
    & WoodLit_coef_water,         &
    & WoodLit_coef_ethanol,       &
    & WoodLit_coef_nonsoluble,    &
    & WoodLit_coef_humus,         &
    & WoodLitterSize,             &
    ! Variable inputs
    & temp,                       &
    & precip,                     &
    & Cflx_wood_2_litter,         &
    & Cflx_green_2_litter,        &
    & Cflx_2_root_exudates,       &
    ! Yasso pools. Used on input and updated on output
    & c_acid_ag1,                 &
    & c_water_ag1,                &
    & c_ethanol_ag1,              &
    & c_nonsoluble_ag1,           &
    & c_acid_bg1,                 &
    & c_water_bg1,                &
    & c_ethanol_bg1,              &
    & c_nonsoluble_bg1,           &
    & c_humus_1,                  &
    & c_acid_ag2,                 &
    & c_water_ag2,                &
    & c_ethanol_ag2,              &
    & c_nonsoluble_ag2,           &
    & c_acid_bg2,                 &
    & c_water_bg2,                &
    & c_ethanol_bg2,              &
    & c_nonsoluble_bg2,           &
    & c_humus_2,                  &
    ! Output fluxes
    & c_decomp_humus_1_sum,       &
    & c_decomp_humus_2_sum,       &
    & c_into_humus_1_sum,         &
    & c_into_humus_2_sum,         &
    & Cflx_soil_respiration       &
    & )

    !$ACC ROUTINE SEQ
    ! Arguments

    REAL(lp), INTENT(in)    :: delta_time_C               !< Carbon time step [s]
    !
    ! Fixed input parameters for Yasso
    !   Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(lp), INTENT(in)    :: LeafLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(lp), INTENT(in)    :: LeafLit_coef_water         ! fraction going to the water soluble pools
    REAL(lp), INTENT(in)    :: LeafLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(lp), INTENT(in)    :: LeafLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(lp), INTENT(in)    :: LeafLit_coef_humus         ! fraction going to the humus pool
    !   Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(lp), INTENT(in)    :: WoodLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(lp), INTENT(in)    :: WoodLit_coef_water         ! fraction going to the water soluble pools
    REAL(lp), INTENT(in)    :: WoodLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(lp), INTENT(in)    :: WoodLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(lp), INTENT(in)    :: WoodLit_coef_humus         ! fraction going to the humus pool
    REAL(lp), INTENT(in)    :: WoodLitterSize             ! size of coarse debris

    ! Carbon inputs [mol m-2 s-1]:
    REAL(lp), INTENT(inout) :: Cflx_wood_2_litter
    REAL(lp), INTENT(inout) :: Cflx_green_2_litter
    REAL(lp), INTENT(inout) :: Cflx_2_root_exudates

    ! Input meteorology for Yasso
    REAL(lp), INTENT(in)    :: temp              ! multi-day mean temperature [K]
    REAL(lp), INTENT(in)    :: precip            ! multi-mean precipitation [mm s-1]

    ! Yasso pools [mol(C) m-2(canopy)], input and updated on output
    ! Size class 1; green litter
    !   Aboveground
    REAL(lp), INTENT(inout) :: c_acid_ag1        ! Acid soluble green litter pool
    REAL(lp), INTENT(inout) :: c_water_ag1       ! Water soluble green litter pool
    REAL(lp), INTENT(inout) :: c_ethanol_ag1     ! Ethanol soluble green litter pool
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1  ! Non-soluble green litter pool
    !   Belowground
    REAL(lp), INTENT(inout) :: c_acid_bg1        ! Acid soluble green litter pool
    REAL(lp), INTENT(inout) :: c_water_bg1       ! Water soluble green litter pool
    REAL(lp), INTENT(inout) :: c_ethanol_bg1     ! Ethanol soluble green litter pool
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1  ! Non-soluble green litter pool
    REAL(lp), INTENT(inout) :: c_humus_1         ! Humus pool from green litter
    ! Size class 2; woody litter
    !   Aboveground
    REAL(lp), INTENT(inout) :: c_acid_ag2        ! Acid soluble woody litter pool
    REAL(lp), INTENT(inout) :: c_water_ag2       ! Water soluble woody litter pool
    REAL(lp), INTENT(inout) :: c_ethanol_ag2     ! Ethanol soluble woody litter pool
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2  ! Non-soluble woody litter pool
    !   Belowground
    REAL(lp), INTENT(inout) :: c_acid_bg2        ! Acid soluble woody litter pool
    REAL(lp), INTENT(inout) :: c_water_bg2       ! Water soluble woody litter pool
    REAL(lp), INTENT(inout) :: c_ethanol_bg2     ! Ethanol soluble woody litter pool
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2  ! Non-soluble woody litter pool
    REAL(lp), INTENT(inout) :: c_humus_2         ! Humus pool from woody litter

    ! Decomposistion fluxes
    REAL(lp), INTENT(inout) :: c_decomp_humus_1_sum   ! annual cflux sum of humus decomposed (leaf)
    REAL(lp), INTENT(inout) :: c_decomp_humus_2_sum   ! annual cflux sum of humus decomposed (wood)
    REAL(lp), INTENT(inout) :: c_into_humus_1_sum     ! annual cflux sum into humus (leaf)
    REAL(lp), INTENT(inout) :: c_into_humus_2_sum     ! annual cflux sum into humus (wood)

    REAL(lp), INTENT(out)   :: Cflx_soil_respiration ! Soil (=heterotrophic) respiration rate  [mol(C) m-2 s-1]
                                                     ! Note: this is a loss, therefore negative!
    !REAL(lp), INTENT(out)   :: soilResp_rate_pot     ! soil respiration without N limitation

    ! local variables for Yasso
    REAL(lp), DIMENSION(18) :: Yasso_io_pools
    REAL(lp), DIMENSION(2)  :: Weather
    REAL(lp), DIMENSION(18) :: Yasso_out
    REAL(lp), DIMENSION(5)  :: LeafLit_coefV
    REAL(lp), DIMENSION(5)  :: WoodLit_coefV
    REAL(lp)                :: leafLitter      ! amount of leafLitter
    REAL(lp)                :: sizeZero

    !REAL(lp) :: Nflux_NPP_2_green_pot,Nflux_NPP_2_wood_pot
    !REAL(lp) :: Nflux_litterGreen_ag_2_slow_pot, Nflux_litterGreen_bg_2_slow_pot
    !REAL(lp) :: Nflux_litter_wood_ag_2_slow_pot,Nflux_litter_wood_bg_2_slow_pot
    !REAL(lp) :: minNflux_litter_green_ag_pot,minNflux_litter_green_bg_pot
    !REAL(lp) :: minNflux_litter_wood_ag_pot,minNflux_litter_wood_bg_pot
    !REAL(lp) :: N2O_emissions_litter_pot

    !REAL(lp) :: minN_plant_demand,minN_soil_demand,total_minN_demand
    !REAL(lp) :: NPP_2_green_via_mobileN, NPP_2_wood_via_mobileN
    !REAL(lp) :: Nmobile_2_green, Nmobile_2_wood, Nmobile_2_green_pot, Nmobile_2_wood_pot
    !REAL(lp) :: redFact_Nmobile
    REAL(lp) :: redFact_Nlimitation
    !REAL(lp) :: Nflx_faeces_2_SMINN, minNflux_litter_active
    !REAL(lp) :: Nflx_crop_harvest_2_SMINN

    !LOGICAL  :: with_Nitrogen        !! .FALSE.: only Carbon pools are updated; .TRUE.: in addition also Nitrogen pools are updated
    !LOGICAL  :: with_yasso           !! .FALSE.: cbalance model for litter and soil carbon is used; .TRUE.: yasso model is used
    !REAL(lp), PARAMETER :: N2O_ef_grazing = 0._lp  !! N2O grazing switched off; reasonable value for the emission factor: 0.0125_lp

    !REAL(lp), PARAMETER :: N2O_ef_grazing = 0._lp  !! N2O grazing switched off; reasonable value for the emission factor: 0.0125_lp

    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    Cflx_soil_respiration     = 0.0_lp
    !soilResp_rate_pot        = 0.0_lp

    ! Preparations
    !-------------------------------------------

    ! TODO: define Yasso_io_pools only with dimension length 9 and set the values anew for the second yasso
    !       call. That way we can use the yasso pool indices defined in mo_carbon_constants.
    !       Perhaps those indices could than be renamed, omitting "POOL_".

    ! Yasso pools:
    Yasso_io_pools(1)  = c_acid_ag1
    Yasso_io_pools(2)  = c_water_ag1
    Yasso_io_pools(3)  = c_ethanol_ag1
    Yasso_io_pools(4)  = c_nonsoluble_ag1
    Yasso_io_pools(5)  = c_acid_bg1
    Yasso_io_pools(6)  = c_water_bg1
    Yasso_io_pools(7)  = c_ethanol_bg1
    Yasso_io_pools(8)  = c_nonsoluble_bg1
    Yasso_io_pools(9)  = c_humus_1
    Yasso_io_pools(10) = c_acid_ag2
    Yasso_io_pools(11) = c_water_ag2
    Yasso_io_pools(12) = c_ethanol_ag2
    Yasso_io_pools(13) = c_nonsoluble_ag2
    Yasso_io_pools(14) = c_acid_bg2
    Yasso_io_pools(15) = c_water_bg2
    Yasso_io_pools(16) = c_ethanol_bg2
    Yasso_io_pools(17) = c_nonsoluble_bg2
    Yasso_io_pools(18) = c_humus_2

    Weather(1) = temp
    Weather(2) = precip

    LeafLit_coefV(LCT_ACID)       = LeafLit_coef_acid
    LeafLit_coefV(LCT_WATER)      = LeafLit_coef_water
    LeafLit_coefV(LCT_ETHANOL)    = LeafLit_coef_ethanol
    LeafLit_coefV(LCT_NONSOLUBLE) = LeafLit_coef_nonsoluble
    LeafLit_coefV(LCT_HUMUS)      = LeafLit_coef_humus
    WoodLit_coefV(LCT_ACID)       = WoodLit_coef_acid
    WoodLit_coefV(LCT_WATER)      = WoodLit_coef_water
    WoodLit_coefV(LCT_ETHANOL)    = WoodLit_coef_ethanol
    WoodLit_coefV(LCT_NONSOLUBLE) = WoodLit_coef_nonsoluble
    WoodLit_coefV(LCT_HUMUS)      = WoodLit_coef_humus

    Yasso_out(1:18)     = 0.0_lp           ! array to store the Yasso pools/fluxes
    sizeZero            = 0.0_lp           ! Green litter by definition has size 0
    redFact_Nlimitation = 1.0_lp           ! Presently no nutrient limitation is used

    !! ============== COMPUTE SOIL RESPIRATION ===================================================================
    Cflx_soil_respiration = 0._lp


    !! call yasso for leaf litter
    CALL yasso(delta_time_C, Yasso_io_pools(1:9), Weather, Cflx_green_2_litter * delta_time_C, &
      &        LeafLit_coefV, sizeZero, Yasso_out, fract_green_aboveGround, &
      &        Cflx_2_root_exudates * delta_time_C, redFact_Nlimitation)

    c_acid_ag1        = Yasso_out(1)
    c_water_ag1       = Yasso_out(2)
    c_ethanol_ag1     = Yasso_out(3)
    c_nonsoluble_ag1  = Yasso_out(4)
    c_acid_bg1        = Yasso_out(5)
    c_water_bg1       = Yasso_out(6)
    c_ethanol_bg1     = Yasso_out(7)
    c_nonsoluble_bg1  = Yasso_out(8)
    c_humus_1         = Yasso_out(9)

    c_decomp_humus_1_sum = c_decomp_humus_1_sum + Yasso_out(17)
    c_into_humus_1_sum   = c_into_humus_1_sum + Yasso_out(18)

    ! rescale respiration (from leaf litter decomposition) from "per day" to "per second"
    Cflx_soil_respiration = Yasso_out(10) / delta_time_C

    !! call yasso for woody litter
    CALL yasso(delta_time_C, Yasso_io_pools(10:18), Weather, Cflx_wood_2_litter * delta_time_C, &
      & WoodLit_coefV, WoodLitterSize, Yasso_out, fract_wood_aboveGround, 0.0_lp, redFact_Nlimitation) ! no exudates

    c_acid_ag2        = Yasso_out(1)
    c_water_ag2       = Yasso_out(2)
    c_ethanol_ag2     = Yasso_out(3)
    c_nonsoluble_ag2  = Yasso_out(4)
    c_acid_bg2        = Yasso_out(5)
    c_water_bg2       = Yasso_out(6)
    c_ethanol_bg2     = Yasso_out(7)
    c_nonsoluble_bg2  = Yasso_out(8)
    c_humus_2         = Yasso_out(9)

    !! Yasso total respiration
    ! add respiration from wood litter decomposition and rescale respiration from "per day" to "per second"
    Cflx_soil_respiration = CFlx_soil_respiration + Yasso_out(10) / delta_time_C

    c_decomp_humus_2_sum = c_decomp_humus_2_sum + Yasso_out(17)
    c_into_humus_2_sum   = c_into_humus_2_sum + Yasso_out(18)

  END SUBROUTINE calc_soil_c


  SUBROUTINE calc_vertical_soil_c(            &
       ! General inputs
    &  delta_time_C, nsoil,                   &
    &  layer_depth, layer_thick,              &
       ! Input coefficients (fixed)
    &  LeafLit_coefV, WoodLit_coefV,          &
    &  MiscPar,                               &
       ! Climate inputs
    &  temp,                                  &
    &  precip,                                &
    &  t_soil_sl,                             &
    &  wtr_soil_rel_sl,                       &
    &  n_thaw_cycle,                          &
       ! Topmodel inputs
    &  frac_inun,                             &
       ! Soil inputs
    &  soil_depth,                            &
    &  soil_depth_sl,                         &
    &  max_depth,                             &
    &  max_depth_sl,                          &
    &  thaw_depth,                            &
    &  permafrost_period_sl,                  &
       ! Litter inputs
    &  woodLitter,                            &
    &  leafLitter,                            &
    &  rootExudates,                          &
       ! Input/output
    &  root_dist_sl,                          &
       ! Yasso pools. Used on input and updated on output
    &  c_acid_ag1,                            &
    &  c_water_ag1,                           &
    &  c_ethanol_ag1,                         &
    &  c_nonsoluble_ag1,                      &
    &  c_acid_ag2,                            &
    &  c_water_ag2,                           &
    &  c_ethanol_ag2,                         &
    &  c_nonsoluble_ag2,                      &
    &  c_acid_bg1_sl,                         &
    &  c_water_bg1_sl,                        &
    &  c_ethanol_bg1_sl,                      &
    &  c_nonsoluble_bg1_sl,                   &
    &  c_humus_1_sl,                          &
    &  c_acid_bg2_sl,                         &
    &  c_water_bg2_sl,                        &
    &  c_ethanol_bg2_sl,                      &
    &  c_nonsoluble_bg2_sl,                   &
    &  c_humus_2_sl,                          &
       ! Output fluxes
    &  c_decomp_humus_1_sum,                  &
    &  c_decomp_humus_2_sum,                  &
    &  c_into_humus_1_sum,                    &
    &  c_into_humus_2_sum,                    &
    &  soilResp_oxic_sl,                      &
    &  soilResp_anox_sl                       &
    &  )

    USE mo_carbon_constants,       ONLY: fract_green_aboveGround, fract_wood_aboveGround,    &
      &                                  LCT_HUMUS, F_RESP, F_2_HUMUS, F_F_HUMUS,            &
      &                                  M_LITTER, M_HUMUS,                                  &
      &                                  WOODLITSIZE, ROOT_DIST_BETA,                        &
      &                                  d_oxic, d_anox, d_size_leaf

    REAL(lp),  INTENT(in)    :: delta_time_C             !< Carbon time step [s]
    INTEGER,   INTENT(in)    :: nsoil                    !
    REAL(lp),  INTENT(in)    :: layer_depth(:)           !< Depths of soil layer center [m]
    REAL(lp),  INTENT(in)    :: layer_thick(:)           !< Thickness of soil layer [m]
    !
    ! Fixed input parameters for Yasso
    !   Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(lp),  INTENT(in)    :: LeafLit_coefV(:)
    !   Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(lp),  INTENT(in)    :: WoodLit_coefV(:)
    !   Other constants
    REAL(lp),  INTENT(in)    :: MiscPar(:)

    ! Carbon inputs [mol m-2 s-1]:
    REAL(lp),  INTENT(in)    :: leafLitter
    REAL(lp),  INTENT(in)    :: woodLitter
    REAL(lp),  INTENT(in)    :: rootExudates

    ! Input meteorology for Yasso
    REAL(lp),  INTENT(in)    :: temp                     !< multi-day mean temperature [K]
    REAL(lp),  INTENT(in)    :: precip                   !< multi-day mean precipitation [mm s-1]
    REAL(lp),  INTENT(in)    :: t_soil_sl(:)             !< soil layer temperatures
    REAL(lp),  INTENT(in)    :: wtr_soil_rel_sl(:)       !< relative soil moisture
    REAL(lp),  INTENT(in)    :: n_thaw_cycle             !< number of freeze-thaw cycles in previous year

    ! Input parameters from topmodel (or other wetland scheme)
    REAL(lp),  INTENT(in)    :: frac_inun                !< Fraction of tile inundated []

    ! Input soil and / or active layer depth
    REAL(lp),  INTENT(in)    :: soil_depth               !< Total soil depth until bedrock [m]
    REAL(lp),  INTENT(in)    :: soil_depth_sl(:)         !< Soil depth within each layer [m]
    REAL(lp),  INTENT(in)    :: max_depth                !< maximum available depth [m], either soil or thaw
    REAL(lp),  INTENT(in)    :: max_depth_sl(:)          !< maximum available depth [m], either soil or thaw, per layer
    REAL(lp),  INTENT(in)    :: thaw_depth               !< thaw depth [m]
    REAL(lp),  INTENT(in)    :: permafrost_period_sl(:)  !< Length of consecutive frozen ground for each soil layer [yr]

    REAL(lp),  INTENT(inout) :: root_dist_sl(:)          !< root distribution, per layer

    ! Yasso pools [mol(C) m-2(canopy)], input and updated on output
    ! Size class 1; green litter
    !   Aboveground
    REAL(lp),  INTENT(inout) :: c_acid_ag1               !< Acid soluble green litter pool
    REAL(lp),  INTENT(inout) :: c_water_ag1              !< Water soluble green litter pool
    REAL(lp),  INTENT(inout) :: c_ethanol_ag1            !< Ethanol soluble green litter pool
    REAL(lp),  INTENT(inout) :: c_nonsoluble_ag1         !< Non-soluble green litter pool
    !   Belowground
    REAL(lp),  INTENT(inout) :: c_acid_bg1_sl(:)         !< Acid soluble green litter pool
    REAL(lp),  INTENT(inout) :: c_water_bg1_sl(:)        !< Water soluble green litter pool
    REAL(lp),  INTENT(inout) :: c_ethanol_bg1_sl(:)      !< Ethanol soluble green litter pool
    REAL(lp),  INTENT(inout) :: c_nonsoluble_bg1_sl(:)   !< Non-soluble green litter pool
    REAL(lp),  INTENT(inout) :: c_humus_1_sl(:)          !< Humus pool from green litter
    ! Size class 2; woody litter
    !   Aboveground
    REAL(lp),  INTENT(inout) :: c_acid_ag2               !< Acid soluble woody litter pool
    REAL(lp),  INTENT(inout) :: c_water_ag2              !< Water soluble woody litter pool
    REAL(lp),  INTENT(inout) :: c_ethanol_ag2            !< Ethanol soluble woody litter pool
    REAL(lp),  INTENT(inout) :: c_nonsoluble_ag2         !< Non-soluble woody litter pool
    !   Belowground
    REAL(lp),  INTENT(inout) :: c_acid_bg2_sl(:)         !< Acid soluble woody litter pool
    REAL(lp),  INTENT(inout) :: c_water_bg2_sl(:)        !< Water soluble woody litter pool
    REAL(lp),  INTENT(inout) :: c_ethanol_bg2_sl(:)      !< Ethanol soluble woody litter pool
    REAL(lp),  INTENT(inout) :: c_nonsoluble_bg2_sl(:)   !< Non-soluble woody litter pool
    REAL(lp),  INTENT(inout) :: c_humus_2_sl(:)          !< Humus pool from woody litter

    ! Decomposition fluxes
    REAL(lp),  INTENT(inout) :: c_decomp_humus_1_sum     !< annual cflux sum of humus decomposed (leaf)
    REAL(lp),  INTENT(inout) :: c_decomp_humus_2_sum     !< annual cflux sum of humus decomposed (wood)
    REAL(lp),  INTENT(inout) :: c_into_humus_1_sum       !< annual cflux sum into humus (leaf)
    REAL(lp),  INTENT(inout) :: c_into_humus_2_sum       !< annual cflux sum into humus (wood)

    REAL(lp),  INTENT(out)   :: soilResp_oxic_sl(:)      !< soil respiration rate, oxic part
    REAL(lp),  INTENT(out)   :: soilResp_anox_sl(:)      !< soil respiration rate, anoxic part

    ! decomposition modifier terms
    REAL(lp)     :: d_temp                  !< term which accounts for the temperature influence on decomposition
    REAL(lp)     :: d_water                 !< term which accounts for the water influence on decomposition
    REAL(lp)     :: d_size_wood             !< term which accounts for the litter size influence on decomposition (woody litter)
    REAL(lp)     :: d_temp_sl(nsoil)        !< term which accounts for the temperature influence on decomposition with depth
    REAL(lp)     :: d_water_sl(nsoil)       !< term which accounts for the water influence on decomposition with depth
    REAL(lp)     :: d_depth_sl(nsoil)       !< term which accounts for the depth influence on decomposition


    REAL(lp) :: frac_inun_loc                   !< local frac_inun -- can be changed without problems

    ! transfer to YASSO subroutine
    REAL(lp), DIMENSION(5)  :: Y_in_pools       !< Yasso pools IN [mol(c)/m2]
    REAL(lp), DIMENSION(2)  :: Y_in_mod         !< Yasso rate modifiers IN
    REAL(lp), DIMENSION(5)  :: Y_out_pools      !< Yasso pools OUT [mol(c)/m2]
    REAL(lp), DIMENSION(4)  :: Y_out_fluxes     !< Yasso fluxes OUT [mol(c)/m2]

    REAL(lp) :: diff_coeff(nsoil)

    REAL(lp) :: tmp_out_humus
    INTEGER  :: is

    REAL(lp) :: soilResp_oxic                   !< ag respiration flux oxic
    REAL(lp) :: soilResp_anox                   !< ag respiration flux anoxic
    CHARACTER(len=*), PARAMETER  :: routine = modname//':calc_vertical_soil_c'


    ! ------------------------------------------------------------
    soilResp_oxic_sl(:)   = 0._lp
    soilResp_anox_sl(:)   = 0._lp

    root_dist_sl = determine_root_dist(max_depth_sl, layer_depth, MiscPar(ROOT_DIST_BETA), t_soil_sl, permafrost_period_sl)

    frac_inun_loc = frac_inun

    ! ------------------------------------------------------------
    ! Formerly within YASSO

    ! ----------------------------------------------------------------------------------------
    ! YASSO calls
    ! ----------------------------------------------------------------------------------------

    ! setup: initialise variables
    tmp_out_humus       = 0._lp
    soilResp_oxic       = 0._lp
    soilResp_anox       = 0._lp
    soilResp_oxic_sl(:) = 0._lp
    soilResp_anox_sl(:) = 0._lp

    ! get decomposition modifiers
    d_temp      = get_modifier_temp(temp)
    d_water     = get_modifier_precip(precip)
    d_temp_sl   = get_modifier_temp(t_soil_sl)
    d_water_sl  = get_modifier_water(wtr_soil_rel_sl)
    d_depth_sl  = get_modifier_depth(layer_depth)
    d_size_wood = get_modifier_size(MiscPar(WOODLITSIZE))

    ! -------------------------------
    ! call YASSO for dry / oxic areas
    !

    ! LEAF LITTER
    ! ag pools

    ! rate modifiers YASSO call for oxic
    Y_in_mod(M_LITTER) = d_temp * d_water * d_size_leaf * d_oxic
    Y_in_mod(M_HUMUS)  = 0._lp   ! no humus above ground

    ! C pools input
    CALL input_YASSO_pools(c_acid_ag1, c_water_ag1, c_ethanol_ag1, &
      &                    c_nonsoluble_ag1, 0._lp, Y_in_pools)
    ! call YASSO for oxic
    CALL yasso_minimal(delta_time_C, Y_in_pools * (1._lp - frac_inun_loc), Y_in_mod, &
      & Y_out_pools=Y_out_pools, Y_out_fluxes=Y_out_fluxes)
    ! C pools output
    CALL output_YASSO_pools(.false., Y_out_pools, c_acid_ag1, c_water_ag1, c_ethanol_ag1, &
      &                     c_nonsoluble_ag1, tmp_out_humus)
    ! Fluxes
    soilResp_oxic              = soilResp_oxic      + Y_out_fluxes(F_RESP) / delta_time_C
    c_into_humus_1_sum         = c_into_humus_1_sum + Y_out_fluxes(F_2_HUMUS)          ! decomposition into humus

    ! call YASSO for inundated / anoxic
    IF (frac_inun_loc > EPSILON(1._lp)) THEN
      Y_in_mod(M_LITTER) = d_temp * d_water * d_size_leaf * d_anox

      ! call YASSO for anoxic
      CALL yasso_minimal(delta_time_C, Y_in_pools * frac_inun_loc, Y_in_mod, &
        & Y_out_pools=Y_out_pools, Y_out_fluxes=Y_out_fluxes)
      ! C pools output
      CALL output_YASSO_pools(.true., Y_out_pools, c_acid_ag1, c_water_ag1, c_ethanol_ag1, &
        &                     c_nonsoluble_ag1, tmp_out_humus)
      ! Fluxes
      soilResp_anox            = soilResp_anox + Y_out_fluxes(F_RESP) / delta_time_C
      c_into_humus_1_sum       = c_into_humus_1_sum   + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus
    END IF

    ! -----------------
    ! YASSO bg1 calls: All soil layers, leaf litter

    DO is = 1, nsoil
      ! rate modifiers YASSO call for oxic
      Y_in_mod(M_LITTER) = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_size_leaf * d_oxic
      Y_in_mod(M_HUMUS)  = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_oxic

      ! C pools input
      CALL input_YASSO_pools(c_acid_bg1_sl(is), c_water_bg1_sl(is), c_ethanol_bg1_sl(is),  &
        &                    c_nonsoluble_bg1_sl(is), c_humus_1_sl(is), Y_in_pools)
      ! call YASSO for oxic
      CALL yasso_minimal(delta_time_C, Y_in_pools * (1._lp - frac_inun_loc), Y_in_mod, &
        & Y_out_pools=Y_out_pools, Y_out_fluxes=Y_out_fluxes)
      ! C pools output
      CALL output_YASSO_pools(.false., Y_out_pools, c_acid_bg1_sl(is), c_water_bg1_sl(is), &
        &                     c_ethanol_bg1_sl(is), c_nonsoluble_bg1_sl(is), c_humus_1_sl(is))
      ! Fluxes
      soilResp_oxic_sl(is)     = soilResp_oxic_sl(is)  + Y_out_fluxes(F_RESP) / delta_time_C     ! respiration
      c_decomp_humus_1_sum     = c_decomp_humus_1_sum  + Y_out_fluxes(F_F_HUMUS)        ! decomposition of humus
      c_into_humus_1_sum       = c_into_humus_1_sum    + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus

      ! call YASSO for inundated / anoxic
      IF (frac_inun_loc > EPSILON(1._lp)) THEN
        Y_in_mod(M_LITTER) = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_size_leaf * d_anox
        Y_in_mod(M_HUMUS)  = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_anox

        ! call YASSO for anoxic
        CALL yasso_minimal(delta_time_C, Y_in_pools * frac_inun_loc, Y_in_mod, &
          & Y_out_pools=Y_out_pools, Y_out_fluxes=Y_out_fluxes)
        ! C pools output
        CALL output_YASSO_pools(.true., Y_out_pools, c_acid_bg1_sl(is), c_water_bg1_sl(is), &
          &                     c_ethanol_bg1_sl(is), c_nonsoluble_bg1_sl(is), c_humus_1_sl(is))
        ! Fluxes
        soilResp_anox_sl(is)   = soilResp_anox_sl(is) + Y_out_fluxes(F_RESP) / delta_time_C      ! respiration
        c_decomp_humus_1_sum   = c_decomp_humus_1_sum + Y_out_fluxes(F_F_HUMUS)        ! decomposition of humus
        c_into_humus_1_sum     = c_into_humus_1_sum   + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus
      END IF
    END DO
    ! some C goes from ag to humus -- accounted for here
    c_humus_1_sl(1) = c_humus_1_sl(1) + tmp_out_humus
    tmp_out_humus   = 0._lp

    ! -----------------
    ! WOODY LITTER
    ! ag pools

    ! rate modifiers YASSO call for oxic
    Y_in_mod(M_LITTER) = d_temp * d_water * d_size_wood * d_oxic
    Y_in_mod(M_HUMUS) = 0._lp   ! no humus above ground

    ! C pools input
    CALL input_YASSO_pools(c_acid_ag2, c_water_ag2, c_ethanol_ag2,  &
      &                    c_nonsoluble_ag2, 0._lp, Y_in_pools)
    ! call YASSO for oxic
    CALL yasso_minimal(delta_time_C, Y_in_pools * (1._lp - frac_inun_loc), Y_in_mod, &
      & Y_out_pools=Y_out_pools, Y_out_fluxes=Y_out_fluxes)
    ! C pools output
    CALL output_YASSO_pools(.false., Y_out_pools, c_acid_ag2, c_water_ag2, c_ethanol_ag2, &
      &                     c_nonsoluble_ag2, tmp_out_humus)
    ! Fluxes
    soilResp_oxic              = soilResp_oxic + Y_out_fluxes(F_RESP) / delta_time_C
    c_into_humus_2_sum         = c_into_humus_2_sum + Y_out_fluxes(F_2_HUMUS)          ! decomposition into humus

    IF (frac_inun_loc > EPSILON(1._lp)) THEN
      ! call YASSO for inundated / anoxic
      Y_in_mod(M_LITTER) = d_temp * d_water * d_size_wood * d_anox

      ! call YASSO for anoxic
      CALL yasso_minimal(delta_time_C, Y_in_pools * frac_inun_loc, Y_in_mod, &
        & Y_out_pools=Y_out_pools, Y_out_fluxes=Y_out_fluxes)
      ! C pools output
      CALL output_YASSO_pools(.true., Y_out_pools, c_acid_ag2, c_water_ag2, c_ethanol_ag2, &
        &                     c_nonsoluble_ag2, tmp_out_humus)
      ! Fluxes
      soilResp_anox            = soilResp_anox + Y_out_fluxes(F_RESP) / delta_time_C
      c_into_humus_2_sum       = c_into_humus_2_sum + Y_out_fluxes(F_2_HUMUS)          ! decomposition into humus
    END IF

    ! -----------------
    ! YASSO bg call: All layers
    DO is = 1, nsoil
      ! rate modifiers, YASSO call for oxic
      Y_in_mod(M_LITTER) = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_size_wood * d_oxic
      Y_in_mod(M_HUMUS)  = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_oxic

      ! C pools input
      CALL input_YASSO_pools(c_acid_bg2_sl(is), c_water_bg2_sl(is), c_ethanol_bg2_sl(is), &
        &                    c_nonsoluble_bg2_sl(is), c_humus_2_sl(is), Y_in_pools)
      ! call YASSO for oxic
      CALL yasso_minimal(delta_time_C, Y_in_pools * (1._lp - frac_inun_loc), Y_in_mod, &
        & Y_out_pools=Y_out_pools, Y_out_fluxes=Y_out_fluxes)
      ! C pools output
      CALL output_YASSO_pools(.false., Y_out_pools, c_acid_bg2_sl(is), c_water_bg2_sl(is), &
        &                     c_ethanol_bg2_sl(is), c_nonsoluble_bg2_sl(is), c_humus_2_sl(is))
      ! Fluxes
      soilResp_oxic_sl(is)     = soilResp_oxic_sl(is) + Y_out_fluxes(F_RESP) / delta_time_C      ! respiration
      c_decomp_humus_2_sum     = c_decomp_humus_2_sum + Y_out_fluxes(F_F_HUMUS)        ! decomposition of humus
      c_into_humus_2_sum       = c_into_humus_2_sum   + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus

      IF (frac_inun_loc > EPSILON(1._lp)) THEN
        ! call YASSO for inundated / anoxic
        Y_in_mod(M_LITTER) = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_size_wood * d_anox
        Y_in_mod(M_HUMUS)  = d_temp_sl(is) * d_depth_sl(is) * d_water_sl(is) * d_anox

        ! call YASSO for anoxic
        CALL yasso_minimal(delta_time_C, Y_in_pools * frac_inun_loc, Y_in_mod, &
          & Y_out_pools=Y_out_pools, Y_out_fluxes=Y_out_fluxes)

        ! C pools output
        CALL output_YASSO_pools(.true., Y_out_pools, c_acid_bg2_sl(is), c_water_bg2_sl(is), &
          &                     c_ethanol_bg2_sl(is), c_nonsoluble_bg2_sl(is), c_humus_2_sl(is))

        ! Fluxes
        soilResp_anox_sl(is)   = soilResp_anox_sl(is) + Y_out_fluxes(F_RESP) / delta_time_C      ! respiration
        c_decomp_humus_2_sum   = c_decomp_humus_2_sum + Y_out_fluxes(F_F_HUMUS)        ! decomposition of humus
        c_into_humus_2_sum     = c_into_humus_2_sum   + Y_out_fluxes(F_2_HUMUS)        ! decomposition into humus
      END IF
    END DO
    ! some C goes from ag to humus -- accounted for here
    c_humus_2_sl(1) = c_humus_2_sl(1) + tmp_out_humus
    tmp_out_humus   = 0._lp

    c_into_humus_1_sum   = c_into_humus_1_sum + &
      &                    LeafLit_coefV(LCT_HUMUS) * leafLitter * delta_time_C
    c_into_humus_2_sum   = c_into_humus_2_sum + &
      &                    WoodLit_coefV(LCT_HUMUS) * woodLitter * delta_time_C

    ! soil respiration: ag respiration is really part of top soil layer respiration
    soilResp_oxic_sl(1)     = soilResp_oxic_sl(1) + soilResp_oxic
    soilResp_anox_sl(1)     = soilResp_anox_sl(1) + soilResp_anox

    ! ----------------------------------------------------------------------------------------
    ! 3b. distribute litter and root exudates
    ! a) leaf litter
    ! ag
    CALL distribute_yasso_litter(c_acid_ag1, c_water_ag1, c_ethanol_ag1, c_nonsoluble_ag1, c_humus_1_sl(1),  &
      &                          leafLitter*delta_time_C, fract_green_aboveGround, LeafLit_coefV)

    ! bg
    DO is = 1, nsoil
      CALL distribute_yasso_litter(c_acid_bg1_sl(is), c_water_bg1_sl(is), c_ethanol_bg1_sl(is),              &
        &                          c_nonsoluble_bg1_sl(is), c_humus_1_sl(is),                                &
        &                          leafLitter*delta_time_C, (1.0_lp - fract_green_aboveGround) * root_dist_sl(is),        &
        &                          LeafLit_coefV)
      ! root exudates are carbohydrates, go directly into water_bg pools
      c_water_bg1_sl(is) = c_water_bg1_sl(is) + rootExudates * delta_time_C * root_dist_sl(is)
    END DO

    !
    ! b) woody litter
    ! ag
    CALL distribute_yasso_litter(c_acid_ag2, c_water_ag2, c_ethanol_ag2, c_nonsoluble_ag2, c_humus_2_sl(1),  &
      &                          woodLitter*delta_time_C, fract_wood_aboveGround, WoodLit_coefV)

    ! bg
    DO is = 1, nsoil
      CALL distribute_yasso_litter(c_acid_bg2_sl(is), c_water_bg2_sl(is), c_ethanol_bg2_sl(is),              &
        &                          c_nonsoluble_bg2_sl(is), c_humus_2_sl(is),                                &
        &                          woodLitter*delta_time_C, (1.0_lp - fract_wood_aboveGround) * root_dist_sl(is),         &
        &                          WoodLit_coefV)
    END DO

    ! ----------------------------------------------------------------------------------------
    ! vertical mixing

    CALL move_ag_2_bg(max_depth, layer_thick, permafrost_period_sl,                                 &
      &               c_acid_ag1, c_acid_bg1_sl, c_acid_ag2, c_acid_bg2_sl,                         &
      &               c_water_ag1, c_water_bg1_sl, c_water_ag2, c_water_bg2_sl,                     &
      &               c_ethanol_ag1, c_ethanol_bg1_sl, c_ethanol_ag2, c_ethanol_bg2_sl,             &
      &               c_nonsoluble_ag1, c_nonsoluble_bg1_sl, c_nonsoluble_ag2, c_nonsoluble_bg2_sl, &
      &               c_humus_1_sl, c_humus_2_sl)

    CALL set_diffusivity(soil_depth, soil_depth_sl, thaw_depth, layer_depth, t_soil_sl, &
      &                  wtr_soil_rel_sl, n_thaw_cycle, permafrost_period_sl, diff_coeff)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_acid_bg1_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_acid_bg2_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_water_bg1_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_water_bg2_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_ethanol_bg1_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_ethanol_bg2_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_nonsoluble_bg1_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_nonsoluble_bg2_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_humus_1_sl)
    CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, c_humus_2_sl)

    ! remove any C that might be below bedrock boundary
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_acid_bg1_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_acid_bg2_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_water_bg1_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_water_bg2_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_ethanol_bg1_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_ethanol_bg2_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_nonsoluble_bg1_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_nonsoluble_bg2_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_humus_1_sl)
    CALL remove_C_below_bedrock(soil_depth, layer_thick, c_humus_2_sl)

  END SUBROUTINE calc_vertical_soil_c

#ifndef _OPENACC
  PURE &
#endif
    SUBROUTINE yasso (delta_time_C, Yasso_in_pools, Weather, c_2_litter, Lit_coefV, WoodLitterSize,       &
    &                 Yasso_out, fract_aboveground, c_2_rootExudates, redFact_Nlimit, is_anoxic, is_peat, &
    &                 A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H)
    !$ACC ROUTINE SEQ
    !! DSG: 11.03.2013
    !! modified Yasso soil C model using an Euler scheme and readable coding style
    !! AWEN (Acid, Water, Ethanol, Nonsoluble) pools are the yasso pools without the humus pools.
    !! JSBACH specific modification of model structure: AWEN pools are separated into aboveground and belowground pools
    !!
    !! Herbivory flux is treated like normal litter, because the herbivory fluxes in JSBACH are not evaluated, yet.
    !! If the represenation of herbivory in JSBACH is evaluated I recommend to treat herbivory flux in YASSO different from general
    !! litter.
    !!
    !! The routine needs as input 1. Yasso pools
    !!                            2. Weather: precipitation and Temperature
    !!                            3. daily litter flux and NPP_2_rootExudates
    !!                            4. PFT specific parameters
    !!
    !! Even when this subroutine is called every day, the dimensions of the variables refere to years!
    !! output pools: AWEN + humus
    !! The routine must be called twice, first for non-woody litter, then for woody litter
    !!-------------------------------------------------------------------------------------------------------------------------
    !! Yasso - Soil carbon model (Jari Liski)
    !!
    !! The code for the yasso model has been developed and is provided by the Finnish Environment Institute SYKE. It is
    !! distributed as part of the ICON Earth System Model under the license conditions as stated in the header of this module.
    !!
    !! model to calculate the amount of soil organic carbon, changes in the amount of soil organic carbon and
    !! heterotrophic soil respiration.
    !! For documention of yasso see Liski et al 2005, Tuomi et al 2008,2009,2011
    !!
    !! implementation: Tea Thum (FMI), Petri Raisanen (FMI), Daniel Goll (MPI-M)
    !!-------------------------------------------------------------------------------------------------------------------------

    USE mo_carbon_constants, ONLY: sec_per_year, LCT_HUMUS,   &
      & POOL_ACID_AG, POOL_WATER_AG, POOL_ETHANOL_AG,         &
      & POOL_NONSOLUBLE_AG,                                   &
      & POOL_ACID_BG, POOL_WATER_BG, POOL_ETHANOL_BG,         &
      & POOL_NONSOLUBLE_BG, POOL_HUMUS,                       &
      & FLUX_RESP_TOTAL, FLUX_2_HUMUS_AG,                     &
      & FLUX_2_HUMUS_BG, FLUX_FROM_HUMUS_DIAG,                &
      & FLUX_RESP_LITTER_AG, FLUX_RESP_LITTER_BG,             &
      & REL_FLUX_AWEN, FLUX_FROM_HUMUS, FLUX_LITTER_2_HUMUS
    USE mo_jsb_physical_constants, ONLY: tmelt
    USE mo_carbon_constants, ONLY: ref_decomp_rate_acid, ref_decomp_rate_water,      &
      & ref_decomp_rate_ethanol, ref_decomp_rate_nonsoluble,  ref_decomp_rate_humus, &
      & temp_p1, temp_p2, precip_p1, size_p1, size_p2, size_p3

    IMPLICIT NONE

    REAL(lp),                INTENT(in)  :: delta_time_C       !< Carbon time step [s]
    ! TODO: instead of Yasso_in_pools and Yasso_out use a variable Yasso_pools or such with intent INOUT
    !       for the Yasso pools and a variable Yasso_fluxes or such replacing Yasso_out(10:18) with intent out.

    REAL(lp), DIMENSION(9),  INTENT(in)  :: Yasso_in_pools     !< Yasso pools in [mol(C) m-2]
    REAL(lp),                INTENT(in)  :: WoodLitterSize     !< Size of woody litter;
                                                               !<   is zero when the subroutine is called for non-woody litter
    REAL(lp), DIMENSION(2),  INTENT(in)  :: Weather            !< climatic drivers: air temperature and precipitation, 15-d means
                                                               !< Weather(1) in [K] ; Weather(2) in [mm s-1]
    REAL(lp),                INTENT(in)  :: c_2_litter         !< fresh litter inputs  (above and belowground) [mol(C) m-2 d-1]
    REAL(lp), DIMENSION(5),  INTENT(in)  :: Lit_coefV          !< fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(lp),                INTENT(in)  :: fract_aboveground  !< parameter to seperate litter influx into above- and belowground
                                                               !<     part                                                   [ ]
    REAL(lp),                INTENT(in)  :: c_2_rootExudates   !< root exudates flux [mol(C) m-2 d-1]

    REAL(lp), DIMENSION(18), INTENT(out) :: Yasso_out          !< updated Y pools (1:9) and fluxes (10:18)
                                                               !<   see update_Cpools for more information

    REAL(lp), OPTIONAL,      INTENT(in)  :: redFact_Nlimit
    LOGICAL,  OPTIONAL,      INTENT(in)  :: is_anoxic
    LOGICAL,  OPTIONAL,      INTENT(in)  :: is_peat

    ! Possibilities to overwrite (some of) the default exchange coefficients (AWEN2H for backward compatibility)
    REAL(lp), OPTIONAL,      INTENT(in)  :: A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H

    ! Fractions to spit the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools,
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(lp)   :: A_2_W
    REAL(lp)   :: A_2_E
    REAL(lp)   :: A_2_N
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(lp)   :: W_2_A
    REAL(lp)   :: W_2_E
    REAL(lp)   :: W_2_N
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(lp)   :: E_2_A
    REAL(lp)   :: E_2_W
    REAL(lp)   :: E_2_N
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(lp)   :: N_2_A
    REAL(lp)   :: N_2_W
    REAL(lp)   :: N_2_E
    ! Alternative pathway for N pool which enters the humus pool [ ]
    REAL(lp)   :: N_2_H
    REAL(lp)   :: AWE_2_H

    ! Yasso pools [mol(C) m-2]
    !   Aboveground
    REAL(lp)     :: c_acid_ag
    REAL(lp)     :: c_water_ag
    REAL(lp)     :: c_ethanol_ag
    REAL(lp)     :: c_nonsoluble_ag
    !   Belowground
    REAL(lp)     :: c_acid_bg
    REAL(lp)     :: c_water_bg
    REAL(lp)     :: c_ethanol_bg
    REAL(lp)     :: c_nonsoluble_bg
    REAL(lp)     :: c_humus

    ! Carbon fluxes [mol(C) m-2 yr-1] - Unit indicated by names CAflx (Carbon annual flux)
    !   Aboveground
    REAL(lp)     :: Cflx_yr_from_acid_ag              !< Loss flux of acid soluble pool
    REAL(lp)     :: Cflx_yr_from_water_ag             !< Loss flux of water soluble pool
    REAL(lp)     :: Cflx_yr_from_ethanol_ag           !< Loss flux of ethanol soluble pool
    REAL(lp)     :: Cflx_yr_from_nonsoluble_ag        !< Loss flux of non soluble pool

    REAL(lp)     :: Cflx_yr_2_acid_ag                 !< Gain flux of acid soluble pool
    REAL(lp)     :: Cflx_yr_2_water_ag                !< Gain flux of water soluble pool
    REAL(lp)     :: Cflx_yr_2_ethanol_ag              !< Gain flux of ethanol soluble pool
    REAL(lp)     :: Cflx_yr_2_nonsoluble_ag           !< Gain flux of non soluble pool

    !   Belowground
    REAL(lp)     :: Cflx_yr_from_acid_bg              !< Loss flux of acid soluble pool
    REAL(lp)     :: Cflx_yr_from_water_bg             !< Loss flux of water soluble pool
    REAL(lp)     :: Cflx_yr_from_ethanol_bg           !< Loss flux of ethanol soluble pool
    REAL(lp)     :: Cflx_yr_from_nonsoluble_bg        !< Loss flux of non soluble pool
    REAL(lp)     :: Cflx_yr_from_humus                !< Loss flux of humus pool

    REAL(lp)     :: Cflx_yr_2_acid_bg                 !< Gain flux of acid soluble pool
    REAL(lp)     :: Cflx_yr_2_water_bg                !< Gain flux of water soluble pool
    REAL(lp)     :: Cflx_yr_2_ethanol_bg              !< Gain flux of ethanol soluble pool
    REAL(lp)     :: Cflx_yr_2_nonsoluble_bg           !< Gain flux of non soluble pool

    !   Humus
    REAL(lp)     :: Cflx_yr_2_humus                   !< Gain flux of humus pool
    REAL(lp)     :: Cflx_yr_2_humusAG                 !< Gain flux of humus pool; from ag only
    REAL(lp)     :: Cflx_yr_2_humusBG                 !< Gain flux of humus pool; from bg only
    REAL(lp)     :: Cflx_yr_litter_2_humus            !< Litter influx to humus soluble pool

    ! Respiration is output of Yasso
    REAL(lp)     :: Cflx_yr_respiration               !< Total flux to the atmosphere
    REAL(lp)     :: Cflx_yr_resp_from_litter_ag       !< Flux to the atmosphere; from AWEN ag only
    REAL(lp)     :: Cflx_yr_resp_from_litter_bg       !< Flux to the atmosphere; from AWEN bg only

    ! Decomposition rates [yr-1]
    REAL(lp)     :: d_acid
    REAL(lp)     :: d_water
    REAL(lp)     :: d_ethanol
    REAL(lp)     :: d_nonsoluble
    REAL(lp)     :: d_humus
    !   For methane relevant processes
    REAL(lp)     :: d_anoxic_ag             !< decomposition rate modifier for ag pools, if conditions are anoxic
    REAL(lp)     :: d_anoxic_bg             !< decomposition rate modifier for bg pools, if conditions are anoxic
    REAL(lp)     :: d_peat_acrotelm         !< decomposition rate modifier for peatland C pools acrotelm
    REAL(lp)     :: d_peat_catotelm         !< decomposition rate modifier for peatland C pools catotelm

    REAL(lp)     :: d_temp                  !< term which accounts for the temperature influence on litter decomposition
    REAL(lp)     :: d_precip                !< term which accounts for the precipitation influence on litter decomposition
    REAL(lp)     :: d_size                  !< term which accounts for the litter size influence on litter decomposition

    ! Decomposition rate for N litter green (needed for N cycle)
    REAL(lp)     :: d_litter_green          !< the rate with which the sum of AWEN pools decays
    REAL(lp)     :: Pseudo_litter_green     !< Sum of all AWEN pools

    ! climatic drivers
    REAL(lp)     :: precip                  !< runnning mean of precipitation      [m yr-1]
    REAL(lp)     :: temp                    !< runnning mean of 2m air temperature [C]

    REAL(lp)     :: redFactor               !< dummy for redFact_Nlimit

    ! anoxic & peat decomposition modifier
    LOGICAL      :: anoxic_flag, peat_flag

    ! time stepping: yasso is parameterized based on annual time
    ! steps, thus the fluxes have to be scaled up to annual rates.
    REAL(lp)     :: years_per_ts, ts_per_year

    ts_per_year = REAL(sec_per_year,lp) / delta_time_C
    years_per_ts = 1.0_lp/ts_per_year

    ! Treatment of optional parameters
    !   Set exchange coefficients
    IF (PRESENT(AWEN2H)) THEN
      AWE_2_H = AWEN2H
      N_2_H   = AWEN2H
    END IF
    ! Default values in [yr-1]
    AWE_2_H  = 0.0045_lp; IF (PRESENT(AWE2H))  AWE_2_H  = AWE2H
    N_2_H    = 0.0045_lp; IF (PRESENT(  N2H))    N_2_H  =   N2H
    A_2_W    = 0.99_lp;   IF (PRESENT(  A2W))    A_2_W  =   A2W
    A_2_E    = 0.00_lp;   IF (PRESENT(  A2E))    A_2_E  =   A2E
    A_2_N    = 0.00_lp;   IF (PRESENT(  A2N))    A_2_N  =   A2N
    W_2_A    = 0.48_lp;   IF (PRESENT(  W2A))    W_2_A  =   W2A
    W_2_E    = 0.00_lp;   IF (PRESENT(  W2E))    W_2_E  =   W2E
    W_2_N    = 0.015_lp;  IF (PRESENT(  W2N))    W_2_N  =   W2N
    E_2_A    = 0.01_lp;   IF (PRESENT(  E2A))    E_2_A  =   E2A
    E_2_W    = 0.00_lp;   IF (PRESENT(  E2W))    E_2_W  =   E2W
    E_2_N    = 0.95_lp;   IF (PRESENT(  E2N))    E_2_N  =   E2N
    N_2_A    = 0.83_lp;   IF (PRESENT(  N2A))    N_2_A  =   N2A
    N_2_W    = 0.01_lp;   IF (PRESENT(  N2W))    N_2_W  =   N2W
    N_2_E    = 0.02_lp;   IF (PRESENT(  N2E))    N_2_E  =   N2E

    anoxic_flag = .FALSE.; IF (PRESENT(is_anoxic)) anoxic_flag = is_anoxic
    peat_flag   = .FALSE.; IF (PRESENT(is_peat))   peat_flag   = is_peat
    redFactor   =  1.0_lp
    IF (PRESENT(redFact_Nlimit)) THEN
      ! nitrogen cycle active & we know already the Nlimitation factor
      redFactor = redFact_Nlimit
    END IF

    ! set up anoxic decomposition rate modifiers
    d_anoxic_ag = 1._lp
    d_anoxic_bg = 1._lp
    IF (anoxic_flag) THEN
      ! anoxic conditions are present, reduce decomposition rates of AWEN pools
      ! following Wania et al. 2009, https://doi.org/10.1029/2008GB003413)
      d_anoxic_ag = 0.35_lp
      d_anoxic_bg = 0.35_lp
    END IF

    ! set up peat decomposition rate modifiers
    d_peat_acrotelm = 1._lp
    d_peat_catotelm = 1._lp
    IF (peat_flag) THEN
      d_peat_acrotelm = 1._lp                   ! no decomposition in litter pools modifier -- we reduce GPP instead
      ! decomposition timescale in LPJ is 7.e-5 -> modifier is applied to (ref_decomp_rate_humus * d_anoxic_bg)
      d_peat_catotelm = 0.125_lp                ! decomposition in catotelm is reduced  -> 0.0016*0.35*0.125 = 7e-5
    END IF

    ! 0. Preparations

    !   Change units of climatic forcing variables
    precip = Weather(2)*sec_per_year/1000._lp  ! mm/s -> m/a
    temp   = Weather(1) - Tmelt                ! K -> C

    !   initialize yasso pools (AWEN + humus)
    c_acid_ag        = Yasso_in_pools(POOL_ACID_AG)
    c_water_ag       = Yasso_in_pools(POOL_WATER_AG)
    c_ethanol_ag     = Yasso_in_pools(POOL_ETHANOL_AG)
    c_nonsoluble_ag  = Yasso_in_pools(POOL_NONSOLUBLE_AG)

    c_acid_bg        = Yasso_in_pools(POOL_ACID_BG)
    c_water_bg       = Yasso_in_pools(POOL_WATER_BG)
    c_ethanol_bg     = Yasso_in_pools(POOL_ETHANOL_BG)
    c_nonsoluble_bg  = Yasso_in_pools(POOL_NONSOLUBLE_BG)

    c_humus          = Yasso_in_pools(POOL_HUMUS)

    !   store sum of the AWEN pools before decomposition
    Pseudo_litter_green = c_acid_ag &
      & + c_water_ag                &
      & + c_ethanol_ag              &
      & + c_nonsoluble_ag           &
      & + c_acid_bg                 &
      & + c_water_bg                &
      & + c_ethanol_bg              &
      & + c_nonsoluble_bg


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 1. Calculate decomposition rates

    !   Litter size dependence of decomposition -- no effect if WoodlitterSize = 0.0
    d_size   = MIN(1.0_lp,(1.0_lp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_lp)**size_p3)
    !   Temperature dependence of decomposition
    d_temp   = EXP(temp_p1*temp + temp_p2*temp**2.0_lp)
    !   Precipitation dependence of decomposition (not for peat)
    IF (peat_flag) THEN
      d_precip = 1.0_lp
    ELSE
      d_precip = 1.0_lp - EXP(precip_p1*precip)
    END IF

    !   decomposition rates accounting for temperature, precipitation, litter size, and nutrient limitation
    d_acid       =  redFactor * ref_decomp_rate_acid       * d_temp * d_precip * d_size * d_peat_acrotelm
    d_water      =  redFactor * ref_decomp_rate_water      * d_temp * d_precip * d_size * d_peat_acrotelm
    d_ethanol    =  redFactor * ref_decomp_rate_ethanol    * d_temp * d_precip * d_size * d_peat_acrotelm
    d_nonsoluble =  redFactor * ref_decomp_rate_nonsoluble * d_temp * d_precip * d_size * d_peat_acrotelm

    !   the decomposition of humus is not limited by nutrients; thus no redFactor
    d_humus      =              ref_decomp_rate_humus      * d_temp * d_precip * d_peat_catotelm ! no size effect on humus

    ! 2. Calculate fluxes

    !   annual loss fluxes (negative values):
    Cflx_yr_from_acid_ag         = c_acid_ag       * d_acid       * d_anoxic_ag
    Cflx_yr_from_water_ag        = c_water_ag      * d_water      * d_anoxic_ag
    Cflx_yr_from_ethanol_ag      = c_ethanol_ag    * d_ethanol    * d_anoxic_ag
    Cflx_yr_from_nonsoluble_ag   = c_nonsoluble_ag * d_nonsoluble * d_anoxic_ag

    Cflx_yr_from_acid_bg         = c_acid_bg       * d_acid       * d_anoxic_bg
    Cflx_yr_from_water_bg        = c_water_bg      * d_water      * d_anoxic_bg
    Cflx_yr_from_ethanol_bg      = c_ethanol_bg    * d_ethanol    * d_anoxic_bg
    Cflx_yr_from_nonsoluble_bg   = c_nonsoluble_bg * d_nonsoluble * d_anoxic_bg

    Cflx_yr_from_humus           = c_humus         * d_humus      * d_anoxic_bg

    !   annual gain fluxes (positive values):
    !   fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why)
    Cflx_yr_2_acid_ag         = ABS(Cflx_yr_from_water_ag      * W_2_A  &    ! returns positive fluxes
      &                           + Cflx_yr_from_ethanol_ag    * E_2_A  &
      &                           + Cflx_yr_from_nonsoluble_ag * N_2_A)

    Cflx_yr_2_water_ag        = ABS(Cflx_yr_from_acid_ag       * A_2_W  &
      &                           + Cflx_yr_from_ethanol_ag    * E_2_W  &
      &                           + Cflx_yr_from_nonsoluble_ag * N_2_W)

    Cflx_yr_2_ethanol_ag      = ABS(Cflx_yr_from_acid_ag       * A_2_E  &
      &                           + Cflx_yr_from_water_ag      * W_2_E  &
      &                           + Cflx_yr_from_nonsoluble_ag * N_2_E)

    Cflx_yr_2_nonsoluble_ag   = ABS(Cflx_yr_from_acid_ag       * A_2_N  &
      &                           + Cflx_yr_from_water_ag      * W_2_N  &
      &                           + Cflx_yr_from_ethanol_ag    * E_2_N)

    Cflx_yr_2_acid_bg         = ABS(Cflx_yr_from_water_bg      * W_2_A  &
      &                           + Cflx_yr_from_ethanol_bg    * E_2_A  &
      &                           + Cflx_yr_from_nonsoluble_bg * N_2_A)

    Cflx_yr_2_water_bg        = ABS(Cflx_yr_from_acid_bg       * A_2_W  &
      &                           + Cflx_yr_from_ethanol_bg    * E_2_W  &
      &                           + Cflx_yr_from_nonsoluble_bg * N_2_W)

    Cflx_yr_2_ethanol_bg      = ABS(Cflx_yr_from_acid_bg       * A_2_E  &
      &                           + Cflx_yr_from_water_bg      * W_2_E  &
      &                           + Cflx_yr_from_nonsoluble_bg * N_2_E)

    Cflx_yr_2_nonsoluble_bg   = ABS(Cflx_yr_from_acid_bg       * A_2_N  &
      &                           + Cflx_yr_from_water_bg      * W_2_N  &
      &                           + Cflx_yr_from_ethanol_bg    * E_2_N)

    Cflx_yr_2_humusAG         = ABS(Cflx_yr_from_acid_ag                &
      &                           + Cflx_yr_from_water_ag               &
      &                           + Cflx_yr_from_ethanol_ag             &
      &                           ) * AWE_2_H                           &
      &                           +                                     &
      &                         ABS(Cflx_yr_from_nonsoluble_ag)         &
      &                           * N_2_H

    Cflx_yr_2_humusBG         = ABS(Cflx_yr_from_acid_bg                &
      &                           + Cflx_yr_from_water_bg               &
      &                           + Cflx_yr_from_ethanol_bg             &
      &                           ) * AWE_2_H                           &
      &                           +                                     &
      &                         ABS(Cflx_yr_from_nonsoluble_bg)         &
      &                           * N_2_H

    Cflx_yr_2_humus           = Cflx_yr_2_humusAG + Cflx_yr_2_humusBG

    !   the remaining fractions of the loss fluxes enter the atmosphere as respiration
!    Cflx_yr_respiration        =  (Cflx_yr_from_acid_ag + Cflx_yr_from_acid_bg)             &
!                              * (1.0_lp - A_2_W - A_2_E - A_2_N - AWE_2_H)            &
!                              + (Cflx_yr_from_water_ag + Cflx_yr_from_water_bg)           &
!                              * (1.0_lp - W_2_A - W_2_E - W_2_N - AWE_2_H)            &
!                              + (Cflx_yr_from_ethanol_ag + Cflx_yr_from_ethanol_bg)       &
!                              * (1.0_lp - E_2_A - E_2_W - E_2_N - AWE_2_H)            &
!                              + (Cflx_yr_from_nonsoluble_ag + Cflx_yr_from_nonsoluble_bg) &
!                              * (1.0_lp - N_2_A - N_2_W - N_2_E - N_2_H)              &
!                              + Cflx_yr_from_humus

    !   litter ag & bg respiration (needed for N cycle)
    Cflx_yr_resp_from_litter_ag = Cflx_yr_from_acid_ag       * (1.0_lp - A_2_W - A_2_E - A_2_N - AWE_2_H) &
      & + Cflx_yr_from_water_ag      * (1.0_lp - W_2_A - W_2_E - W_2_N - AWE_2_H) &
      & + Cflx_yr_from_ethanol_ag    * (1.0_lp - E_2_A - E_2_W - E_2_N - AWE_2_H) &
      & + Cflx_yr_from_nonsoluble_ag * (1.0_lp - N_2_A - N_2_W - N_2_E -   N_2_H)

    Cflx_yr_resp_from_litter_bg = Cflx_yr_from_acid_bg       * (1.0_lp - A_2_W - A_2_E - A_2_N - AWE_2_H) &
      & + Cflx_yr_from_water_bg      * (1.0_lp - W_2_A - W_2_E - W_2_N - AWE_2_H) &
      & + Cflx_yr_from_ethanol_bg    * (1.0_lp - E_2_A - E_2_W - E_2_N - AWE_2_H) &
      & + Cflx_yr_from_nonsoluble_bg * (1.0_lp - N_2_A - N_2_W - N_2_E -   N_2_H)

    Cflx_yr_respiration         = Cflx_yr_resp_from_litter_ag                     &
      & + Cflx_yr_resp_from_litter_bg                                             &
      & + Cflx_yr_from_humus

    ! 3. update Yasso pools

    !   3a. redistribution among yasso pools
    c_acid_ag        = MAX(0.0_lp,                                         &
      &                    c_acid_ag                                       & ! old pool
      &                    + (Cflx_yr_from_acid_ag                         & ! incoming flux from AWEN pools
      &                    + Cflx_yr_2_acid_ag) * years_per_ts)
    c_water_ag       = MAX(0.0_lp,                                         & ! and so on .....
      &                    c_water_ag                                      &
      &                    + (Cflx_yr_from_water_ag                        &
      &                    + Cflx_yr_2_water_ag) * years_per_ts)
    c_ethanol_ag     = MAX(0.0_lp,                                         &
      &                    c_ethanol_ag                                    &
      &                    + (Cflx_yr_from_ethanol_ag                      &
      &                    + Cflx_yr_2_ethanol_ag) * years_per_ts)
    c_nonsoluble_ag  = MAX(0.0_lp,                                         &
      &                    c_nonsoluble_ag                                 &
      &                    + (Cflx_yr_from_nonsoluble_ag                   &
      &                    + Cflx_yr_2_nonsoluble_ag) * years_per_ts)

    c_acid_bg        = MAX(0.0_lp,                                         &
      &                    c_acid_bg                                       & ! old pool
      &                    + (Cflx_yr_from_acid_bg                         & ! incoming flux from AWEN pools
      &                    + Cflx_yr_2_acid_bg) * years_per_ts)
    c_water_bg       = MAX(0.0_lp,                                         & ! ...
      &                    c_water_bg                                      &
      &                    + (Cflx_yr_from_water_bg                        &
      &                    + Cflx_yr_2_water_bg) * years_per_ts            &
      & + c_2_rootExudates)                               ! exudates are carbohydrates only and belowground
    c_ethanol_bg     = MAX(0.0_lp,                                         &
      &                    c_ethanol_bg                                    &
      &                    + (Cflx_yr_from_ethanol_bg                      &
      &                    + Cflx_yr_2_ethanol_bg) * years_per_ts)
    c_nonsoluble_bg  = MAX(0.0_lp,                                         &
      &                    c_nonsoluble_bg                                 &
      &                    + (Cflx_yr_from_nonsoluble_bg                   &
      &                    + Cflx_yr_2_nonsoluble_bg) * years_per_ts)

    c_humus          = MAX(0.0_lp,                                         &
      &                    c_humus                                         &
      &                    + (Cflx_yr_from_humus                           &
      &                    + Cflx_yr_2_humus) * years_per_ts)

    !   3b. distribute litter
    !    ag
    CALL distribute_yasso_litter(c_acid_ag, c_water_ag, c_ethanol_ag, c_nonsoluble_ag, c_humus, &
      &                          c_2_litter,          fract_aboveground, Lit_coefV)
    !    bg
    CALL distribute_yasso_litter(c_acid_bg, c_water_bg, c_ethanol_bg, c_nonsoluble_bg, c_humus, &
      &                          c_2_litter, 1.0_lp - fract_aboveground, Lit_coefV)

    ! Cflx_litter_2_humus is not only required for distribution, but possibly also as output for the eq script
    Cflx_yr_litter_2_humus = Lit_coefV(LCT_HUMUS) * c_2_litter * ts_per_year

    ! compute d_litter_green; this is given by the change in AWEN pools
    !          Flx                 Pool(t)                   Pool(t+1)
    !    d_P = ----     &   Flx = ----------   --->   d_P = -----------  - 1
    !          Pool(t)             Pool(t+1)                 Pool(t)

    IF (Pseudo_litter_green > 0.0_lp) THEN
       d_litter_green = (  c_acid_ag          &
         &               + c_water_ag         &
         &               + c_ethanol_ag       &
         &               + c_nonsoluble_ag    &
         &               + c_acid_bg          &
         &               + c_water_bg         &
         &               + c_ethanol_bg       &
         &               + c_nonsoluble_bg)   &
         & / Pseudo_litter_green - 1._lp
    ELSE
      d_litter_green = 0.0_lp
    END IF

    ! 4.1 Write pools into the output array
    Yasso_out(POOL_ACID_AG)       = c_acid_ag
    Yasso_out(POOL_WATER_AG)      = c_water_ag
    Yasso_out(POOL_ETHANOL_AG)    = c_ethanol_ag
    Yasso_out(POOL_NONSOLUBLE_AG) = c_nonsoluble_ag
    Yasso_out(POOL_ACID_BG)       = c_acid_bg
    Yasso_out(POOL_WATER_BG)      = c_water_bg
    Yasso_out(POOL_ETHANOL_BG)    = c_ethanol_bg
    Yasso_out(POOL_NONSOLUBLE_BG) = c_nonsoluble_bg
    Yasso_out(POOL_HUMUS)         = c_humus

    ! 4.2 Write respiration into the output array
    Yasso_out(FLUX_RESP_TOTAL)    = Cflx_yr_respiration * years_per_ts ! convert back to daily

    IF (.NOT. PRESENT(redFact_Nlimit)) THEN
      ! the routine might have been called to diagnose the nutrient demand, so we give out
      ! 1. flux for immobilisation: this is the sum of all carbon coming from the non-humus pools to humus
      Yasso_out(FLUX_2_HUMUS_AG)      =  Cflx_yr_2_humusAG * years_per_ts !  per timestep ! the flux is positive
      Yasso_out(FLUX_2_HUMUS_BG)      =  Cflx_yr_2_humusBG * years_per_ts !  per timestep ! the flux is positive
      ! 2. flux for mineralisation:  this is the flux from humus to atmosphere
      Yasso_out(FLUX_FROM_HUMUS_DIAG) = -Cflx_yr_from_humus * years_per_ts !  per timestep ! sign CHECKED!
      ! ... as well as all fluxes from AWEN pools to atmosphere
      Yasso_out(FLUX_RESP_LITTER_AG)  = -Cflx_yr_resp_from_litter_ag * years_per_ts ! per timestep ! the flux is negative
      Yasso_out(FLUX_RESP_LITTER_BG)  = -Cflx_yr_resp_from_litter_bg * years_per_ts ! per timestep ! the flux is negative
      ! 3. decomposition rate for Npool litter green
      Yasso_out(REL_FLUX_AWEN)        =  d_litter_green
    END IF

    ! 4.3 Write diagnostics to determine equilibrium of humus pools
    Yasso_out(FLUX_FROM_HUMUS)     = -Cflx_yr_from_humus * years_per_ts
    Yasso_out(FLUX_LITTER_2_HUMUS) = (Cflx_yr_2_humus + Cflx_yr_litter_2_humus) * years_per_ts

  END SUBROUTINE yasso

  ! ====================================================================================================== !
  !
  !> Distributes incomming litter to yasso pools (i.e. proportional to litter coefficients and passed factors)
  !
#ifndef _OPENACC
  PURE &
#endif
    SUBROUTINE distribute_yasso_litter(c_acid, c_water, c_ethanol, c_nonsoluble, c_humus, &
    &                                  litter, fraction, coefficients)

    !$ACC ROUTINE SEQ

    USE mo_carbon_constants, ONLY:   LCT_ACID, LCT_WATER, LCT_ETHANOL, LCT_NONSOLUBLE, LCT_HUMUS

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    REAL(lp), INTENT(inout)             :: c_acid        !< yasso acid pool
    REAL(lp), INTENT(inout)             :: c_water       !< yasso water pool
    REAL(lp), INTENT(inout)             :: c_ethanol     !< yasso ethanol pool
    REAL(lp), INTENT(inout)             :: c_nonsoluble  !< yasso nonsoluble pool
    REAL(lp), INTENT(inout)             :: c_humus       !< yasso humus pool
    REAL(lp), INTENT(in)                :: litter        !< litter input (above and belowground) [mol(C) m-2 day-1]
    REAL(lp), INTENT(in)                :: fraction      !< fraction of litter put into yasso pools []
    REAL(lp), DIMENSION(5), INTENT(in)  :: coefficients  !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !

    CALL add_litter_to_yasso_pool(c_acid,       litter, fraction, coefficients(LCT_ACID))
    CALL add_litter_to_yasso_pool(c_water,      litter, fraction, coefficients(LCT_WATER))
    CALL add_litter_to_yasso_pool(c_ethanol,    litter, fraction, coefficients(LCT_ETHANOL))
    CALL add_litter_to_yasso_pool(c_nonsoluble, litter, fraction, coefficients(LCT_NONSOLUBLE))
    CALL add_litter_to_yasso_pool(c_humus,      litter, fraction, coefficients(LCT_HUMUS))

  END SUBROUTINE distribute_yasso_litter

  ! ====================================================================================================== !
  !
  !> Adds incomming litter share (i.e. proportional to litter coefficient and factor)
  !
#ifndef _OPENACC
  ELEMENTAL &
#endif
    SUBROUTINE add_litter_to_yasso_pool(this_pool, litter, fraction, coefficient)
    !$ACC ROUTINE SEQ
    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    REAL(lp), INTENT(inout) :: this_pool     !< target yasso pool
    REAL(lp), INTENT(in)    :: litter        !< litter input  (above and belowground) [mol(C) m-2 day-1]
    REAL(lp), INTENT(in)    :: fraction      !< fraction of litter put into this pool (BG vs AG) [ ]
    REAL(lp), INTENT(in)    :: coefficient   !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !
    this_pool           = MAX(0.0_lp, this_pool + (fraction * coefficient * litter))

  END SUBROUTINE add_litter_to_yasso_pool

  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2
  !
  SUBROUTINE get_per_tile(  &
    & Cpool_per_tile,       & ! out
    & Cpool_per_canopy,     & ! in
    & veg_fract_correction, & ! in
    & fract_fpc_max)          ! in

    REAL(lp), INTENT(in), DIMENSION(:) :: &
      & Cpool_per_canopy, &
      & fract_fpc_max,    &
      & veg_fract_correction
    REAL(lp), INTENT(out), DIMENSION(:) :: &
      & Cpool_per_tile

    INTEGER :: ic

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, SIZE(fract_fpc_max)
      Cpool_per_tile(ic) = Cpool_per_canopy(ic) * veg_fract_correction(ic) * fract_fpc_max(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE get_per_tile
  SUBROUTINE yasso_minimal(delta_time_C, Y_in_pools, Y_in_mod,                                             &
    &                      A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H, &
    &                      Y_out_pools, Y_out_fluxes)

    USE mo_carbon_constants, ONLY: Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS, &
      &  F_RESP, F_2_HUMUS, F_F_HUMUS, F_LITRESP, ref_decomp_rate_acid,               &
      &  ref_decomp_rate_water, ref_decomp_rate_ethanol, ref_decomp_rate_nonsoluble,  &
      &  ref_decomp_rate_humus, M_LITTER, M_HUMUS, sec_per_year

    REAL(lp), INTENT(in)                 :: delta_time_C     !< Carbon time step [s]
    REAL(lp), DIMENSION(5),  INTENT(in)  :: Y_in_pools       !< Yasso pools IN [mol(c)/m2]
    REAL(lp), DIMENSION(2),  INTENT(in)  :: Y_in_mod         !< Yasso rate modifiers IN

    ! Possibilities to overwrite (some of) the default exchange coefficients (AWEN2H for backward compatibility)
    REAL(lp), OPTIONAL,      INTENT(in)  :: A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H
    REAL(lp), DIMENSION(5),  INTENT(out) :: Y_out_pools      !< Yasso pools OUT [mol(c)/m2]
    REAL(lp), DIMENSION(4),  INTENT(out) :: Y_out_fluxes     !< Yasso fluxes OUT [mol(c)/m2]

    !local
    ! Yasso pools   [mol(c)/m2]
    REAL(lp)     :: c_acid
    REAL(lp)     :: c_water
    REAL(lp)     :: c_ethanol
    REAL(lp)     :: c_nonsoluble
    REAL(lp)     :: c_humus

    ! Decomposition rates [1/a]
    REAL(lp)     :: d_acid
    REAL(lp)     :: d_water
    REAL(lp)     :: d_ethanol
    REAL(lp)     :: d_nonsoluble
    REAL(lp)     :: d_humus

    ! Carbon fluxes [mol(C) m-2 yr-1] - Unit indicated by names Cflx_yr (Carbon annual flux)
    REAL(lp)     :: Cflx_yr_from_acid              !< Loss flux of acid soluble pool        [mol(c)/m2/yr]
    REAL(lp)     :: Cflx_yr_from_water             !< Loss flux of water soluble pool       [mol(c)/m2/yr]
    REAL(lp)     :: Cflx_yr_from_ethanol           !< Loss flux of ethanol soluble pool     [mol(c)/m2/yr]
    REAL(lp)     :: Cflx_yr_from_nonsoluble        !< Loss flux of non soluble pool         [mol(c)/m2/yr]
    REAL(lp)     :: Cflx_yr_from_humus             !< Loss flux of humus pool               [mol(c)/m2/yr]

    REAL(lp)     :: Cflx_yr_2_acid                 !< Gain flux of acid soluble pool        [mol(c)/m2/yr]
    REAL(lp)     :: Cflx_yr_2_water                !< Gain flux of water soluble pool       [mol(c)/m2/yr]
    REAL(lp)     :: Cflx_yr_2_ethanol              !< Gain flux of ethanol soluble pool     [mol(c)/m2/yr]
    REAL(lp)     :: Cflx_yr_2_nonsoluble           !< Gain flux of non soluble pool         [mol(c)/m2/yr]
    REAL(lp)     :: Cflx_yr_2_humus                !< Gain flux of humus pool               [mol(c)/m2/yr]

    REAL(lp)     :: Cflx_yr_resp_from_litter
    REAL(lp)     :: Cflx_yr_respiration

    ! Fractions to split the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools,
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(lp)   :: A_2_W
    REAL(lp)   :: A_2_E
    REAL(lp)   :: A_2_N
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(lp)   :: W_2_A
    REAL(lp)   :: W_2_E
    REAL(lp)   :: W_2_N
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(lp)   :: E_2_A
    REAL(lp)   :: E_2_W
    REAL(lp)   :: E_2_N
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(lp)   :: N_2_A
    REAL(lp)   :: N_2_W
    REAL(lp)   :: N_2_E
    ! fraction of decomposition fluxes of AWEN pools which enters the humus pool [ ]
    ! REAL(lp)   :: AWEN_2_H
    ! Alternative pathway for N pool which enters the humus pool [ ]
    REAL(lp)   :: N_2_H
    REAL(lp)   :: AWE_2_H

    ! time stepping; the parameterization of yasso is done on a annual time
    ! step, thus the fluxes have to be scaled up to annual rates.
    REAL(lp)     :: years_per_ts, ts_per_year

    ts_per_year = REAL(sec_per_year,lp) / delta_time_C
    years_per_ts = 1.0_lp/ts_per_year

    ! Treatment of optional parameters
    !   Set exchange coefficients
    IF (PRESENT(AWEN2H)) THEN
      AWE_2_H = AWEN2H
      N_2_H   = AWEN2H
    END IF

    ! Default values in [yr-1]
    AWE_2_H  = 0.0045_lp; IF (PRESENT(AWE2H))  AWE_2_H  = AWE2H
    N_2_H    = 0.0045_lp; IF (PRESENT(  N2H))    N_2_H  =   N2H
    A_2_W    = 0.99_lp;   IF (PRESENT(  A2W))    A_2_W  =   A2W
    A_2_E    = 0.00_lp;   IF (PRESENT(  A2E))    A_2_E  =   A2E
    A_2_N    = 0.00_lp;   IF (PRESENT(  A2N))    A_2_N  =   A2N
    W_2_A    = 0.48_lp;   IF (PRESENT(  W2A))    W_2_A  =   W2A
    W_2_E    = 0.00_lp;   IF (PRESENT(  W2E))    W_2_E  =   W2E
    W_2_N    = 0.015_lp;  IF (PRESENT(  W2N))    W_2_N  =   W2N
    E_2_A    = 0.01_lp;   IF (PRESENT(  E2A))    E_2_A  =   E2A
    E_2_W    = 0.00_lp;   IF (PRESENT(  E2W))    E_2_W  =   E2W
    E_2_N    = 0.95_lp;   IF (PRESENT(  E2N))    E_2_N  =   E2N
    N_2_A    = 0.83_lp;   IF (PRESENT(  N2A))    N_2_A  =   N2A
    N_2_W    = 0.01_lp;   IF (PRESENT(  N2W))    N_2_W  =   N2W
    N_2_E    = 0.02_lp;   IF (PRESENT(  N2E))    N_2_E  =   N2E

    ! INPUTS
    ! a) pools
    c_acid       = Y_in_pools(Y_ACID)
    c_water      = Y_in_pools(Y_WATER)
    c_ethanol    = Y_in_pools(Y_ETHANOL)
    c_nonsoluble = Y_in_pools(Y_NONSOLUBLE)
    c_humus      = Y_in_pools(Y_HUMUS)

    ! b) decomposition modifiers
    ! decomposition rate is reference decomposition rate * modifier
    ! with modifier function(temperature, water, size, etc) set externally
    d_acid       = ref_decomp_rate_acid       * Y_in_mod(M_LITTER)
    d_water      = ref_decomp_rate_water      * Y_in_mod(M_LITTER)
    d_ethanol    = ref_decomp_rate_ethanol    * Y_in_mod(M_LITTER)
    d_nonsoluble = ref_decomp_rate_nonsoluble * Y_in_mod(M_LITTER)
    d_humus      = ref_decomp_rate_humus      * Y_in_mod(M_HUMUS)

    ! Calculate fluxes
    !   annual loss fluxes (negative values):
    Cflx_yr_from_acid         = c_acid       * d_acid
    Cflx_yr_from_water        = c_water      * d_water
    Cflx_yr_from_ethanol      = c_ethanol    * d_ethanol
    Cflx_yr_from_nonsoluble   = c_nonsoluble * d_nonsoluble
    Cflx_yr_from_humus        = c_humus      * d_humus

    !   annual gain fluxes (positive values):
    !   fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why)
    Cflx_yr_2_acid            = ABS(Cflx_yr_from_water      * W_2_A  &    ! returns positive fluxes
      &                           + Cflx_yr_from_ethanol    * E_2_A  &
      &                           + Cflx_yr_from_nonsoluble * N_2_A)

    Cflx_yr_2_water           = ABS(Cflx_yr_from_acid       * A_2_W  &
      &                            + Cflx_yr_from_ethanol    * E_2_W &
      &                            + Cflx_yr_from_nonsoluble * N_2_W)

    Cflx_yr_2_ethanol         = ABS(Cflx_yr_from_acid       * A_2_E  &
      &                            + Cflx_yr_from_water      * W_2_E &
      &                            + Cflx_yr_from_nonsoluble * N_2_E)

    Cflx_yr_2_nonsoluble      = ABS(Cflx_yr_from_acid       * A_2_N  &
      &                           + Cflx_yr_from_water      * W_2_N  &
      &                           + Cflx_yr_from_ethanol    * E_2_N)

    Cflx_yr_2_humus           = ABS(Cflx_yr_from_acid                &
      &                           + Cflx_yr_from_water               &
      &                           + Cflx_yr_from_ethanol             &
      &                           ) * AWE_2_H                        &
      &                             +                                &
      &                         ABS(Cflx_yr_from_nonsoluble)         &
      &                             * N_2_H

    !   litter ag & bg respiration (needed for N cycle)
    Cflx_yr_resp_from_litter =  Cflx_yr_from_acid       * (1.0_lp - A_2_W - A_2_E - A_2_N - AWE_2_H) &
      &                       + Cflx_yr_from_water      * (1.0_lp - W_2_A - W_2_E - W_2_N - AWE_2_H) &
      &                       + Cflx_yr_from_ethanol    * (1.0_lp - E_2_A - E_2_W - E_2_N - AWE_2_H) &
      &                       + Cflx_yr_from_nonsoluble * (1.0_lp - N_2_A - N_2_W - N_2_E -   N_2_H)

    Cflx_yr_respiration      =  Cflx_yr_resp_from_litter                      &
      &                       + Cflx_yr_from_humus

    ! 3. update Yasso pools

    !   3a. redistribution among yasso pools
    c_acid        = MAX(0.0_lp,                                    &
      &                 c_acid                                     & ! old pool
      &              + (Cflx_yr_from_acid                          & ! incoming flux from AWEN pools
      &               + Cflx_yr_2_acid) * years_per_ts)
    c_water       = MAX(0.0_lp,                                    & ! ...
      &                 c_water                                    &
      &              + (Cflx_yr_from_water                         &
      &               + Cflx_yr_2_water) * years_per_ts)
    c_ethanol     = MAX(0.0_lp,                                    &
      &                 c_ethanol                                  &
      &              + (Cflx_yr_from_ethanol                       &
      &               + Cflx_yr_2_ethanol) * years_per_ts)
    c_nonsoluble  = MAX(0.0_lp,                                    &
      &                 c_nonsoluble                               &
      &              + (Cflx_yr_from_nonsoluble                    &
      &               + Cflx_yr_2_nonsoluble) * years_per_ts)

    c_humus       = MAX(0.0_lp,                                    &
      &                 c_humus                                    &
      &              + (Cflx_yr_from_humus                         &
      &               + Cflx_yr_2_humus) * years_per_ts)

    ! 4.1 Write pools into the output array
    Y_out_pools(Y_ACID)       = c_acid
    Y_out_pools(Y_WATER)      = c_water
    Y_out_pools(Y_ETHANOL)    = c_ethanol
    Y_out_pools(Y_NONSOLUBLE) = c_nonsoluble
    Y_out_pools(Y_HUMUS)      = c_humus

    ! 4.2 Write respiration into the output array
    Y_out_fluxes(F_RESP)     = Cflx_yr_respiration * years_per_ts ! convert back to daily

    ! the routine might have been called to diagnose the nutrient demand, so we give out
    ! 1. flux for immobilisation: this is the sum of all carbon coming from the non-humus pools to humus
    Y_out_fluxes(F_2_HUMUS)  =  Cflx_yr_2_humus * years_per_ts !  per timestep ! the flux is positive
    ! 2. flux for mineralisation:  this is the flux from humus to atmosphere
    Y_out_fluxes(F_F_HUMUS)  = -Cflx_yr_from_humus * years_per_ts !  per timestep ! sign CHECKED!
    ! ... as well as all fluxes from AWEN pools to atmosphere
    Y_out_fluxes(F_LITRESP)  = -Cflx_yr_resp_from_litter * years_per_ts ! per timestep ! the flux is negative
  END SUBROUTINE yasso_minimal

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> We here determine the distribution of plant roots in the soil (based on Jackson 1996), assuming an
  !> exponential decay of root material with depth.
  !>
  !> Note: The root distribution assumed here is only used to determine the vertical profile of root litter components.
  !>       It has no direct affect on the plant hydrology.
  !
  FUNCTION determine_root_dist(max_thick, layer_depth, root_dist_beta, t_soil_sl, permafrost_period_sl) RESULT(root_dist)
    USE mo_carbon_constants,       ONLY: permafrost_threshold
    USE mo_jsb_physical_constants, ONLY: tmelt

    REAL(lp), INTENT(in)  :: max_thick(:)            !< Active layer thickness [m]
    REAL(lp), INTENT(in)  :: layer_depth(:)          !< Model layer depth [m]
    REAL(lp), INTENT(in)  :: root_dist_beta          !< PFT dependent factor for root decay with depth
    REAL(lp), INTENT(in)  :: t_soil_sl(:)            !< Soil temperature [K]
    REAL(lp), INTENT(in)  :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]


    ! locals
    CHARACTER(len=*), PARAMETER :: routine = modname//':determine_root_dist'
    REAL(lp) :: root_dist(size(max_thick))     ! root distribution in each layer
    REAL(lp) :: root_dist_cum(size(max_thick)) ! cumulative root distribution for layers
    INTEGER  :: is, nsoil, nsoil_p

    nsoil = size(max_thick)
    nsoil_p = nsoil
    root_dist(:) = 0._lp
    root_dist_cum(:) = 0._lp

    ! Determine which soil layers belong to the active layer depth limited by either permafrost
    ! diagnostics or temperatures below zero.
    DO is = 1, nsoil
      IF ((permafrost_period_sl(is) > permafrost_threshold) .OR. (t_soil_sl(is) < tmelt)) THEN
        nsoil_p = MAX(is-1, 1)
        EXIT
      END IF
    END DO

    IF (root_dist_beta > 0._lp) THEN
      DO is = 1, nsoil_p
        root_dist_cum(is) = 1._lp - (root_dist_beta ** (layer_depth(is) * 100._lp))
          IF (is == 1) THEN
            root_dist(is) = root_dist_cum(is)
          ELSE
            root_dist(is) = root_dist_cum(is) - root_dist_cum(is-1)
          END IF
      END DO
    END IF

    ! No root distribution is computed for permafrost affected soil layers
    DO is = 1, nsoil
      IF (permafrost_period_sl(is) > permafrost_threshold) THEN
        root_dist(is) = 0._lp
      END IF
    END DO

    ! normalise root distribution -> SUM(root_dist) = 1
    IF (SUM(root_dist) > 0._lp) THEN
      root_dist(:) = root_dist(:) / SUM(root_dist(:))
    ELSE
      root_dist(:) = 0._lp
    END IF

  END FUNCTION determine_root_dist

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> The routine prepares the yasso call: The "yasso_in" buffer is filled with data from the current vertical
  !> yasso pools.
  !>
  PURE SUBROUTINE input_YASSO_pools(cpool_acid, cpool_water, cpool_ethanol,      &
    &                               cpool_nonsoluble, cpool_humus, YASSO_in)
    USE mo_carbon_constants,  ONLY: Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS

    REAL(lp), INTENT(in)  :: cpool_acid
    REAL(lp), INTENT(in)  :: cpool_water
    REAL(lp), INTENT(in)  :: cpool_ethanol
    REAL(lp), INTENT(in)  :: cpool_nonsoluble
    REAL(lp), INTENT(in)  :: cpool_humus
    REAL(lp), INTENT(out) :: YASSO_in(5)

    YASSO_in(Y_ACID)       = cpool_acid
    YASSO_in(Y_WATER)      = cpool_water
    YASSO_in(Y_ETHANOL)    = cpool_ethanol
    YASSO_in(Y_NONSOLUBLE) = cpool_nonsoluble
    YASSO_in(Y_HUMUS)      = cpool_humus
  END SUBROUTINE input_YASSO_pools

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> This subroutine moves the yasso output to the respective vertical yasso pools. Depending on argument "plus"
  !> the yasso output values are added to (true) or are replacing (false) current vertical pool values.
  !>
  PURE SUBROUTINE output_YASSO_pools(plus, YASSO_out, cpool_acid, cpool_water, cpool_ethanol, &
    &                                cpool_nonsoluble, cpool_humus)

    USE mo_carbon_constants,  ONLY: Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS

    LOGICAL,  INTENT(in)     :: plus
    REAL(lp), INTENT(in)     :: YASSO_out(5)
    REAL(lp), INTENT(inout)  :: cpool_acid
    REAL(lp), INTENT(inout)  :: cpool_water
    REAL(lp), INTENT(inout)  :: cpool_ethanol
    REAL(lp), INTENT(inout)  :: cpool_nonsoluble
    REAL(lp), INTENT(inout)  :: cpool_humus

    IF (plus) THEN
      cpool_acid       = cpool_acid       + YASSO_out(Y_ACID)
      cpool_water      = cpool_water      + YASSO_out(Y_WATER)
      cpool_ethanol    = cpool_ethanol    + YASSO_out(Y_ETHANOL)
      cpool_nonsoluble = cpool_nonsoluble + YASSO_out(Y_NONSOLUBLE)
      cpool_humus      = cpool_humus      + YASSO_out(Y_HUMUS)
    ELSE
      cpool_acid       = YASSO_out(Y_ACID)
      cpool_water      = YASSO_out(Y_WATER)
      cpool_ethanol    = YASSO_out(Y_ETHANOL)
      cpool_nonsoluble = YASSO_out(Y_NONSOLUBLE)
      cpool_humus      = YASSO_out(Y_HUMUS)
    END IF
  END SUBROUTINE output_YASSO_pools

  ! ---------------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION get_modifier_temp(temp_soil, use_Q10) RESULT(d_temp)
    USE mo_carbon_constants,       ONLY: Q10, temp_p1, temp_p2
    USE mo_jsb_physical_constants, ONLY: tmelt

    REAL(lp), INTENT(in)          :: temp_soil
    LOGICAL, OPTIONAL, INTENT(in) :: use_Q10

    LOGICAL :: do_Q10
    REAL(lp):: d_temp

    d_temp = 0._lp

    do_Q10 = .FALSE.
    IF (PRESENT(use_Q10)) THEN
      do_Q10 = use_Q10
    END IF

    IF (.NOT. do_Q10) THEN       ! use YASSO formulation
      d_temp  = EXP(temp_p1*(temp_soil - tmelt) + &                 ! YASSO formulation for temperature dependence
        &           temp_p2*(temp_soil- tmelt)**2.0_lp)
    ELSE                        ! use alternative Q10 formulation
      d_temp  = Q10**((temp_soil - (tmelt + 10._lp))/10._lp)        ! Temperature modifier as CLM
    END IF
  END FUNCTION get_modifier_temp

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_water(water_soil) RESULT(d_water)
    USE mo_carbon_constants,  ONLY: moist_p_a, moist_p_b, moist_p_d

    REAL(lp), INTENT(in)  :: water_soil(:)

    REAL(lp) :: d_water(size(water_soil))
    INTEGER  :: is, nsoil

    nsoil = size(water_soil)
    d_water(:) = 0._lp

    DO is = 1,nsoil
      d_water(is) = ( ( (water_soil(is) - moist_p_b) / (moist_p_a-moist_p_b) )**  &
        &             ( moist_p_d*( (moist_p_b-moist_p_a)/moist_p_a ) ) ) *       &
        &             ( (water_soil(is) / moist_p_a)**moist_p_d )
    END DO
  END FUNCTION get_modifier_water

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_depth(layer_depth) RESULT(d_depth)
    USE mo_carbon_constants,  ONLY: ref_dep

    REAL(lp), INTENT(in)  :: layer_depth(:)

    REAL(lp) :: d_depth(size(layer_depth))
    INTEGER :: is, nsoil
    nsoil = size(layer_depth)
    d_depth(:) = 0._lp

    DO is = 1, nsoil
      d_depth(is) = exp(-layer_depth(is)/ref_dep)                       ! depth modifier as CLM
    END DO

  END FUNCTION get_modifier_depth

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_size(WoodLitterSize) RESULT(d_size)
    USE mo_carbon_constants,       ONLY: size_p1, size_p2, size_p3

    REAL(lp), INTENT(in)  :: WoodLitterSize
    REAL(lp) :: d_size

    d_size = MIN(1.0_lp,(1.0_lp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_lp)**size_p3)

  END FUNCTION get_modifier_size

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_precip(precip) RESULT(d_precip)
    USE mo_carbon_constants,       ONLY: precip_p1, sec_per_year

    REAL(lp), INTENT(in)  :: precip
    REAL(lp) :: d_precip

    d_precip = 1.0_lp - EXP(precip_p1*(precip*sec_per_year/1000._lp))      ! YASSO formulation for precipitation dependence

  END FUNCTION get_modifier_precip

  ! ---------------------------------------------------------------------------------------------------------
  SUBROUTINE vertical_C_transport(act_depth_sl, diff_coeff, delta_time_C, dz, permafrost_period_sl, cpool_bg_sl)
    USE mo_carbon_constants,       ONLY: minimum_layer_thickness, permafrost_threshold

    REAL(lp), INTENT(in)    :: act_depth_sl(:)         !< actual layer thickness, considering soil and active layer depth [m]
    REAL(lp), INTENT(in)    :: diff_coeff(:)           !< diffusivity
    REAL(lp), INTENT(in)    :: delta_time_C            !< Carbon time step [s]
    REAL(lp), INTENT(in)    :: dz(:)                   !< layer thickness [m]
    REAL(lp), INTENT(in)    :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(inout) :: cpool_bg_sl(:)          !< below ground C pools on soil layers

    ! layer counters, # soil layers
    INTEGER :: is, nsoil, nsoil_loc
    REAL(lp), PARAMETER :: implicit_contrib = 0.5_lp !< weighting factor for Crank-Nicholson scheme

    ! Diffusion constants averaged over adjacent layers.
    REAL(lp) :: D_below, D_above

    ! Layer-dependent weighting factors in Crank-Nicholson scheme.
    REAL(lp) :: dz_factor, W_below, W_above
    REAL(lp) :: dz_here, dz_above, dz_below
    REAL(lp) :: C_here, C_above, C_below

    ! Diagonal, left and right subdiagonals for Crank-Nicholson matrix
    REAL(lp) :: left(size(act_depth_sl))
    REAL(lp) :: diag(size(act_depth_sl))
    REAL(lp) :: right(size(act_depth_sl))

    ! Requires carbon concentrations for diffusion algorithm
    REAL(lp) :: C_conc(size(act_depth_sl))
    REAL(lp) :: src_C_conc(size(act_depth_sl))

    ! Right hand side vector for Crank-Nicholson scheme equations
    ! This is overwritten with solution by dgtsv tridiagonal solver
    REAL(lp) :: rhs(size(act_depth_sl))
    INTEGER  :: info             ! dgtsv status variable

    REAL(lp) :: C_before, C_after

    nsoil = size(act_depth_sl)
    nsoil_loc = nsoil
    src_C_conc(:) = 0._lp

    C_before = SUM(cpool_bg_sl(:))

    IF (C_before > 0._lp) THEN
      C_conc(:) = 0._lp
      DO is = 1, nsoil
        IF ((act_depth_sl(is) > minimum_layer_thickness) .AND. (permafrost_period_sl(is) < permafrost_threshold)) THEN
          ! only do this, if layer is thick enough and non-permafrost
          C_conc(is) = cpool_bg_sl(is) / act_depth_sl(is)
        ELSE
          nsoil_loc = MAX(is - 1, 1)
          EXIT
        END IF
      END DO

      !------------------------------------------------------------
      ! build tridiagonal matrix and known right hand side

      ! Process the local layers
      DO is = 1, nsoil_loc

        ! Extract sensible values to use for concentration and
        ! thickness of the current layer, the layer above and below
        dz_here = dz(is)
        C_here = C_conc(is)
        IF (is == 1) THEN
          dz_above = 0._lp
          C_above = 0._lp
        ELSE
          dz_above = dz(is-1)
          C_above = C_conc(is-1)
        END IF
        IF (is == nsoil_loc) THEN
          dz_below = 0._lp
          C_below = 0._lp
        ELSE
          dz_below = dz(is+1)
          C_below = C_conc(is+1)
        END IF

        ! Determine diffusion coefficients to the layers above and below the current layer.
        IF (is == 1) THEN
          ! upwards diffusivity in uppermost layer is zero
          D_above = 0._lp
        ELSE
          ! for the layers below, upwards diffusivity is the weighted mean of the layer
          ! diffusivity in layer and layer above
          D_above = (dz(is-1)*diff_coeff(is-1) + &
            &        dz(is  )*diff_coeff(is))  / &
            &       (dz(is-1) + dz(is))
        END IF
        IF (is == nsoil_loc) THEN
          ! downwards diffusivity in the bottom layer is set to zero
          ! to impose no flux boundary condition
          D_below = 0._lp
        ELSE
          ! for the layers above, downwards diffusivity is the weighted mean of the
          ! diffusivity in the layer and the layer below
          D_below = (dz(is  )*diff_coeff(is  )  + &
            &        dz(is+1)*diff_coeff(is+1)) / &
            &       (dz(is) + dz(is+1))
        END IF

        ! Calculate layer-dependent weighting factors.
        dz_factor = 0.25_lp * (dz_below + 2._lp * dz_here + dz_above)
        W_below = D_below * (delta_time_C / dz_factor) * (2._lp / (dz_below + dz_here))
        W_above = D_above * (delta_time_C / dz_factor) * (2._lp / (dz_here + dz_above))

        ! boundary conditions for surface and bottom layers: no flux through top / bottom
        IF (is == 1) THEN
          W_above = 0._lp
        END IF
        IF (is == nsoil_loc) THEN
          W_below = 0._lp
        END IF

        ! Fill in matrix diagonal and off-diagonal elements.
        diag(is) = 1.0_lp + (implicit_contrib * (W_below + W_above))
        left(is) = -implicit_contrib * W_above
        right(is) = -implicit_contrib * W_below

        ! Calculate right hand side vector values.
        rhs(is) = (1.0_lp - (1.0_lp - implicit_contrib)*(W_below + W_above)) * C_here + &
          &        (1.0_lp - implicit_contrib)*W_below*C_below +                         &
          &        (1.0_lp - implicit_contrib)*W_above*C_above +                         &
          &        src_C_conc(is)
      END DO

      ! Semi-implicit case -- whenever implicit_contrib > 0.1,
      ! fully explicit case is solved already ( = rhs)
      IF (implicit_contrib >= 0.1_lp) THEN
        !------------------------------------------------------------
        ! solve tridiagonal system using dgtsv from NETLIB
        CALL dgtsv(nsoil_loc, 1, left(2:nsoil_loc), diag, &
          &        right(1:nsoil_loc-1), rhs, nsoil_loc, info)
  !       IF (PRESENT(lissue)) lissue = info /= 0
  !      IF (info /= 0) THEN
  !        WRITE(message_text,*) 'dgtsv info ', info
  !        CALL message(routine, message_text, all_print=.TRUE., level=em_warn)
  !      END IF
      END IF

      ! solution is ALWAYS written into the rhs vector, also in the fully explicit case
      C_conc = rhs

      DO is = 1, nsoil_loc
        IF ((act_depth_sl(is) > minimum_layer_thickness) .AND. (permafrost_period_sl(is) < permafrost_threshold)) THEN
          ! only do this, if layer is thick enough and non-permafrost
          cpool_bg_sl(is) = C_conc(is) * act_depth_sl(is)
        END IF
      END DO

      WHERE(cpool_bg_sl(:) < 0._lp)
        cpool_bg_sl(:) = 0._lp
      END WHERE

      C_after = SUM(cpool_bg_sl(:))

      ! re-normalise to conserve carbon
      IF (C_after > 0._lp) THEN
        cpool_bg_sl(:) = cpool_bg_sl(:)/C_after*C_before
      END IF
    END IF    ! C_before > 0

  END SUBROUTINE vertical_C_transport

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> Computes the diffusivity for the vertical carbon transport through the soil column. Permafrost
  !> regions are dominated by cryoturbation and the diffusivity computation is limited to the
  !> active layer depth. In regions not affected by permafrost, the vertical transport is
  !> facilitated by bioturbation.
  !>
  SUBROUTINE set_diffusivity(soil_depth, soil_depth_sl, thaw_depth, layer_depth, t_soil_sl, &
    &                        wtr_rel_sl, n_thaw_cycle, permafrost_period_sl, diff_coeff)
    USE mo_jsb_physical_constants, ONLY: tmelt
    USE mo_carbon_constants,       ONLY: minimum_layer_thickness, permafrost_threshold

    REAL(lp), INTENT(in)    :: soil_depth
    REAL(lp), INTENT(in)    :: soil_depth_sl(:)
    REAL(lp), INTENT(in)    :: thaw_depth
    REAL(lp), INTENT(in)    :: layer_depth(:)
    REAL(lp), INTENT(in)    :: t_soil_sl(:)
    REAL(lp), INTENT(in)    :: wtr_rel_sl(:)
    REAL(lp), INTENT(in)    :: n_thaw_cycle
    REAL(lp), INTENT(in)    :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(out)   :: diff_coeff(:)
    INTEGER :: is, isp1, nsoil

    REAL(lp), PARAMETER :: cryo_max = 4.7e-12_lp !< Maximum diffusivity for cryoturbation [m^2 s-1]
    REAL(lp), PARAMETER :: bio_max  = 4.7e-13_lp !< Maximum diffusivity for bioturbation [m^2 s-1]


    nsoil = size(diff_coeff)

    diff_coeff(:) = 0._lp

    IF (thaw_depth < soil_depth) THEN    ! active layer depth < soil depth --> permafrost
      DO is = 1,nsoil
        isp1 = MIN(is+1,nsoil)
        IF (  (permafrost_period_sl(is)   < permafrost_threshold) .AND. &
          &   (permafrost_period_sl(isp1) < permafrost_threshold) .AND. &
          &   (t_soil_sl(is)          > tmelt)) THEN
          diff_coeff(is) = cryo_max *                               &
            &             MIN(1._lp, wtr_rel_sl(is)) *              &
            &              MIN(1._lp, (n_thaw_cycle / 40._lp))*     &
            &              (1._lp - layer_depth(is)/thaw_depth)
        ELSE
          diff_coeff(is) = 0._lp
        END IF
      END DO
    ELSE
      DO is = 1,nsoil
        isp1 = MIN(is+1,nsoil)
        IF ((soil_depth_sl(is) > minimum_layer_thickness) .AND. (t_soil_sl(is) > tmelt)) THEN
          diff_coeff(is) = MAX(0._lp, bio_max *(1._lp - (layer_depth(is)/3._lp)))
        ELSE
          diff_coeff(is) = 0._lp
        END IF
      END DO
    END IF

  END SUBROUTINE set_diffusivity

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> This subroutine moves carbon from the above ground litter pools to the below ground litter pools,
  !> if the above ground litter pools exceed a certain threshold. The threshold is defined as a fraction
  !> of the maximum carbon content in the above ground layer. The excess carbon is moved to the uppermost
  !> below ground layer. Similarily, if the below ground litter pools exceed a certain threshold, the
  !> excess carbon is moved to the next below ground layer.
  !>
    SUBROUTINE move_ag_2_bg(max_depth, layer_thick, permafrost_period_sl,                                     &
      &                     c_acid_ag1, c_acid_bg1_sl, c_acid_ag2, c_acid_bg2_sl,                         &
      &                     c_water_ag1, c_water_bg1_sl, c_water_ag2, c_water_bg2_sl,                     &
      &                     c_ethanol_ag1, c_ethanol_bg1_sl, c_ethanol_ag2, c_ethanol_bg2_sl,             &
      &                     c_nonsoluble_ag1, c_nonsoluble_bg1_sl, c_nonsoluble_ag2, c_nonsoluble_bg2_sl, &
      &                     c_humus_1_sl, c_humus_2_sl)

    USE mo_carbon_constants, ONLY: C_2_biomass_ratio, rho_org_ag, rho_org_bg, molarMassC_kg, &
      &                            permafrost_threshold

    REAL(lp), INTENT(in)    :: max_depth               !< Soil depth until bedrock or permafrost [m]
    REAL(lp), INTENT(in)    :: layer_thick(:)          !< Thickness of soil layers [m]
    REAL(lp), INTENT(in)    :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(inout) :: c_acid_ag1              !< Above ground acid soluble green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg1_sl(:)        !< Below ground acid soluble green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_ag2              !< Above ground acid soluble woody litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_acid_bg2_sl(:)        !< Below ground acid soluble woody litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_ag1             !< Above ground water soluble green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg1_sl(:)       !< Below ground water soluble green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_ag2             !< Above ground water soluble woody litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_water_bg2_sl(:)       !< Below ground water soluble woody litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_ag1           !< Above ground ethanol soluble green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1_sl(:)     !< Below ground ethanol soluble green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_ag2           !< Above ground ethanol soluble woody litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2_sl(:)     !< Below ground ethanol soluble woody litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1        !< Above ground nonsoluble green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1_sl(:)  !< Below ground nonsoluble green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2        !< Above ground nonsoluble woody litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2_sl(:)  !< Below ground nonsoluble woody litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_1_sl(:)         !< Humus pool from green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_2_sl(:)         !< Humus pool from woody litter pool [mol(C) m-2]

    REAL(lp),  PARAMETER :: thresh_frac_org    = 1.25_lp  !< Threshold organic fraction above which C is being moved downwards
    REAL(lp),  PARAMETER :: ag_layer_depth     = 0.13_lp  !< Assumed depth of above ground litter layer [m]
    REAL(lp) :: cpool_tot
    REAL(lp) :: depth_upper
    REAL(lp) :: rfrac
    REAL(lp) :: tmp

    REAL(lp) :: max_C_ag
    REAL(lp) :: max_C_bg(size(layer_thick))
    INTEGER  :: is, isp1, nsoil

    nsoil = SIZE(layer_thick)

    ! determine max C amount per layer
    ! C [kgC] = height [m] x A [m2] x rho [kg/m3] x C_2_bio [kgC/kg]
    ! finally converted to mol
    max_C_ag = ag_layer_depth * rho_org_ag * C_2_biomass_ratio / molarMassC_kg
    max_C_bg(:) = layer_thick(:) * rho_org_bg * C_2_biomass_ratio / molarMassC_kg

    ! Check ag soil
    cpool_tot = c_acid_ag1 + c_acid_ag2       +      &
      &         c_water_ag1 + c_water_ag2     +      &
      &         c_ethanol_ag1 + c_ethanol_ag2 +      &
      &         c_nonsoluble_ag1 + c_nonsoluble_ag2

    IF (cpool_tot > max_C_ag) THEN
      ! Rescale above ground carbon pools such, that maximum amount is no longer exceeded
      ! and move excess carbon into the uppermost below ground layer.

      rfrac = max_C_ag / cpool_tot

      tmp                    = c_acid_ag1
      c_acid_ag1             = c_acid_ag1 * rfrac
      c_acid_bg1_sl(1)       = c_acid_bg1_sl(1) + (tmp - c_acid_ag1)
      tmp                    = c_acid_ag2
      c_acid_ag2             = c_acid_ag2 * rfrac
      c_acid_bg2_sl(1)       = c_acid_bg2_sl(1) + (tmp - c_acid_ag2)
      tmp                    = c_water_ag1
      c_water_ag1            = c_water_ag1 * rfrac
      c_water_bg1_sl(1)      = c_water_bg1_sl(1) + (tmp - c_water_ag1)
      tmp                    = c_water_ag2
      c_water_ag2            = c_water_ag2 * rfrac
      c_water_bg2_sl(1)      = c_water_bg2_sl(1) + (tmp - c_water_ag2)
      tmp                    = c_ethanol_ag1
      c_ethanol_ag1          = c_ethanol_ag1 * rfrac
      c_ethanol_bg1_sl(1)    = c_ethanol_bg1_sl(1) + (tmp - c_ethanol_ag1)
      tmp                    = c_ethanol_ag2
      c_ethanol_ag2          = c_ethanol_ag2 * rfrac
      c_ethanol_bg2_sl(1)    = c_ethanol_bg2_sl(1) + (tmp - c_ethanol_ag2)
      tmp                    = c_nonsoluble_ag1
      c_nonsoluble_ag1       = c_nonsoluble_ag1 * rfrac
      c_nonsoluble_bg1_sl(1) = c_nonsoluble_bg1_sl(1) + (tmp - c_nonsoluble_ag1)
      tmp                    = c_nonsoluble_ag2
      c_nonsoluble_ag2       = c_nonsoluble_ag2 * rfrac
      c_nonsoluble_bg2_sl(1) = c_nonsoluble_bg2_sl(1) + (tmp - c_nonsoluble_ag2)
    END IF

    ! check bg soil layers, move C downwards, if possible
    depth_upper = 0._lp
    DO is = 1,nsoil-1
      isp1 = MIN(is+1, nsoil)
      IF ((depth_upper < max_depth) .AND. (permafrost_period_sl(is) < permafrost_threshold) &
        & .AND. (permafrost_period_sl(isp1) < permafrost_threshold)) THEN
        depth_upper = depth_upper + layer_thick(is)
        cpool_tot = c_acid_bg1_sl(is)       + c_acid_bg2_sl(is)       +      &
          &         c_water_bg1_sl(is)      + c_water_bg2_sl(is)      +      &
          &         c_ethanol_bg1_sl(is)    + c_ethanol_bg2_sl(is)    +      &
          &         c_nonsoluble_bg1_sl(is) + c_nonsoluble_bg2_sl(is) +      &
          &         c_humus_1_sl(is)        + c_humus_2_sl(is)

        IF (cpool_tot > max_C_bg(is)) THEN
          ! Rescale this layers carbon pools such, that maximum amount is no longer exceeded
          ! and move excess carbon into the layer below.

          rfrac = max_C_bg(is) / cpool_tot

          tmp                       = c_acid_bg1_sl(is)
          c_acid_bg1_sl(is)         = c_acid_bg1_sl(is) * rfrac
          c_acid_bg1_sl(is+1)       = c_acid_bg1_sl(is+1) + (tmp - c_acid_bg1_sl(is))
          tmp                       = c_acid_bg2_sl(is)
          c_acid_bg2_sl(is)         = c_acid_bg2_sl(is) * rfrac
          c_acid_bg2_sl(is+1)       = c_acid_bg2_sl(is+1) + (tmp - c_acid_bg2_sl(is))
          tmp                       = c_water_bg1_sl(is)
          c_water_bg1_sl(is)        = c_water_bg1_sl(is) * rfrac
          c_water_bg1_sl(is+1)      = c_water_bg1_sl(is+1) + (tmp - c_water_bg1_sl(is))
          tmp                       = c_water_bg2_sl(is)
          c_water_bg2_sl(is)        = c_water_bg2_sl(is) * rfrac
          c_water_bg2_sl(is+1)      = c_water_bg2_sl(is+1) + (tmp - c_water_bg2_sl(is))
          tmp                       = c_ethanol_bg1_sl(is)
          c_ethanol_bg1_sl(is)      = c_ethanol_bg1_sl(is) * rfrac
          c_ethanol_bg1_sl(is+1)    = c_ethanol_bg1_sl(is+1) + (tmp - c_ethanol_bg1_sl(is))
          tmp                       = c_ethanol_bg2_sl(is)
          c_ethanol_bg2_sl(is)      = c_ethanol_bg2_sl(is) * rfrac
          c_ethanol_bg2_sl(is+1)    = c_ethanol_bg2_sl(is+1) + (tmp - c_ethanol_bg2_sl(is))
          tmp                       = c_nonsoluble_bg1_sl(is)
          c_nonsoluble_bg1_sl(is)   = c_nonsoluble_bg1_sl(is) * rfrac
          c_nonsoluble_bg1_sl(is+1) = c_nonsoluble_bg1_sl(is+1) + (tmp - c_nonsoluble_bg1_sl(is))
          tmp                       = c_nonsoluble_bg2_sl(is)
          c_nonsoluble_bg2_sl(is)   = c_nonsoluble_bg2_sl(is) * rfrac
          c_nonsoluble_bg2_sl(is+1) = c_nonsoluble_bg2_sl(is+1) + (tmp - c_nonsoluble_bg2_sl(is))
          tmp                       = c_humus_1_sl(is)
          c_humus_1_sl(is)          = c_humus_1_sl(is) * rfrac
          c_humus_1_sl(is+1)        = c_humus_1_sl(is+1) + (tmp - c_humus_1_sl(is))
          tmp                       = c_humus_2_sl(is)
          c_humus_2_sl(is)          = c_humus_2_sl(is) * rfrac
          c_humus_2_sl(is+1)        = c_humus_2_sl(is+1) + (tmp - c_humus_2_sl(is))
        END IF
      END IF
    END DO
  END SUBROUTINE move_ag_2_bg

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> Remove any carbon that is below the bedrock depth and move it to the lowest layer above bedrock.
  !> This should not happen in reality, but can occur due to numerical issues.
  !> This routine is a safeguard against that.
  !>
  SUBROUTINE remove_C_below_bedrock(soil_depth, layer_thick, c_pool_bg)
    REAL(lp), INTENT(in)     :: soil_depth
    REAL(lp), INTENT(in)     :: layer_thick(:)
    REAL(lp), INTENT(inout)  :: c_pool_bg(:)

    INTEGER  :: is, nsoil, max_layer_above
    REAL(lp) :: soil_above, C_below_bedrock, upper_bound

    nsoil = size(layer_thick)

    soil_above = 0._lp
    C_below_bedrock = 0._lp
    max_layer_above = 0
    upper_bound = 0._lp

    ! Check for C below bedrock depth
    DO is = 1, nsoil
      IF (soil_depth > upper_bound) THEN
        max_layer_above = is
        soil_above = soil_above + layer_thick(is)
      ELSE
        IF (c_pool_bg(is) > 0._lp) THEN
          C_below_bedrock = C_below_bedrock + c_pool_bg(is)
        END IF
      END IF
      upper_bound = upper_bound + layer_thick(is)
    END DO

    ! Move C to lowest layer above bedrock, the lowest layer with a soil fraction
    IF ((C_below_bedrock > 0._lp) .AND. (soil_above > 0._lp) .AND. (max_layer_above > 0)) THEN
      c_pool_bg(max_layer_above) = c_pool_bg(max_layer_above) + C_below_bedrock
      DO is = MIN(max_layer_above+1,nsoil), nsoil
        c_pool_bg(is) = 0._lp
      END DO
    END IF

  END SUBROUTINE remove_C_below_bedrock

  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of fire for the carbon pools. More precisely:
  ! It is assumed that for the burned area the carbon from the above ground litter pools (c_litter_green_ag,
  ! c_litter_wood_ag) is completely released to the atmosphere and the carbon from the living plant pools
  ! (c_green_st, c_reserve_st, c_woods_st) is partly released to the atmosphere and partly relocated
  ! into the litter pools.
  !
  SUBROUTINE relocate_carbon_fire( &
    & delta_time_C,                      & ! in
    & burned_fract,                      & ! in from disturbed_frac
    & fire_fract_wood_2_atmos,           & ! in
    & LeafLit_coef,                      & ! in
    & WoodLit_coef,                      & ! in
    & fuel_moisture,                     & ! in
    & c_green,                           & ! inout
    & c_woods,                           & ! inout
    & c_reserve,                         & ! inout
    & cflux_dist_green_2_soil,           & ! inout
    & cflux_dist_woods_2_soil,           & ! inout
    & c_acid_ag1,                        & ! inout
    & c_water_ag1,                       & ! inout
    & c_ethanol_ag1,                     & ! inout
    & c_nonsoluble_ag1,                  & ! inout
    & c_acid_ag2,                        & ! inout
    & c_water_ag2,                       & ! inout
    & c_ethanol_ag2,                     & ! inout
    & c_nonsoluble_ag2,                  & ! inout
    & c_acid_bg1,                        & ! inout
    & c_water_bg1,                       & ! inout
    & c_ethanol_bg1,                     & ! inout
    & c_nonsoluble_bg1,                  & ! inout
    & c_acid_bg2,                        & ! inout
    & c_water_bg2,                       & ! inout
    & c_ethanol_bg2,                     & ! inout
    & c_nonsoluble_bg2,                  & ! inout
    & c_humus_1,                         & ! inout
    & c_humus_2,                         & ! inout
    & co2flux_fire_all_2_atm             & ! out
    & )

    ! Arguments
    REAL(lp), INTENT(in)    :: delta_time_C                 !< Carbon time step [s]
    REAL(lp), INTENT(in)    :: burned_fract(:)              ! Fraction of the vegetated area of each tile burned till the last
                                                            ! call of this routine
    REAL(lp), INTENT(in)    :: fire_fract_wood_2_atmos      ! Fraction of above ground wood immediately emitted to atmosphere
                                                            ! by fire
    REAL(lp), INTENT(in)    :: fuel_moisture(:)             !
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)              ! Factor to spread non woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)              ! Factor to spread woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(inout) :: c_green(:)                   ! Value of green carbon pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_woods(:)                   ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_reserve(:)                 ! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(out)   :: co2flux_fire_all_2_atm(:)    ! CO2 immediately released by fire per tile [kg(CO2)m-2(canopy)s-1]
    REAL(lp), INTENT(inout) :: cflux_dist_green_2_soil(:)   ! Amount of carbon relocated by wind and fire damage
                                                            ! .. to the green litter pools [mol(C)/m^2(canopy)s-1]
    REAL(lp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   ! Amount of carbon relocated by wind damage and fire
                                                            ! .. to the wood litter pools [mol(C)/m^2(canopy)s-1]
    !   YASSO
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    REAL(lp), INTENT(inout) :: c_acid_ag1(:)          ! Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_ag1(:)         ! Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_ag1(:)       ! Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1(:)    ! Non-soluble green litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp), INTENT(inout) :: c_acid_bg1(:)          ! Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg1(:)         ! Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1(:)       ! Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1(:)    ! Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_1(:)           ! Humus pool from green litter [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    REAL(lp), INTENT(inout) :: c_acid_ag2(:)          ! Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_ag2(:)         ! Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_ag2(:)       ! Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2(:)    ! Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp), INTENT(inout) :: c_acid_bg2(:)          ! Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg2(:)         ! Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2(:)       ! Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2(:)    ! Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_2(:)           ! Humus pool from woody litter [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(lp)  :: non_woody_litter, fraction ! Helper
    INTEGER   :: ic, nc

    nc = SIZE(burned_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc
      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Amount of carbon released from the green and reserve pools entering the green litter pools
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + &
        & (c_green(ic) + c_reserve(ic)) * (1._lp - fract_green_aboveGround) * burned_fract(ic) / delta_time_C

      ! Amount of carbon released from the wood pool entering the woody litter pools
      cflux_dist_woods_2_soil(ic)  = cflux_dist_woods_2_soil(ic) + &
        & c_woods(ic) * (1._lp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * burned_fract(ic) / delta_time_C

      ! Determin the amount of carbon released from the wood pool, the living tissue pools, and the
      ! litter pools going to the atmosphere
      co2flux_fire_all_2_atm(ic) = ( &
          ! above ground fraction of living biomass
        & non_woody_litter * fract_green_aboveGround                          &
          ! add the above ground YASSO pools
        & + c_water_ag1(ic) + c_acid_ag1(ic) + c_ethanol_ag1(ic)              &
        & + c_nonsoluble_ag1(ic)                                              &
        & + c_water_ag2(ic) + c_acid_ag2(ic) + c_ethanol_ag2(ic)              &
        & + c_nonsoluble_ag2(ic)                                              &
          ! add the fire affected fraction of the above ground wood pools
        & + c_woods(ic) * fire_fract_wood_2_atmos * fract_wood_aboveGround    &
          ! scale with the burned fraction
        & ) * burned_fract(ic) * molarMassCO2_kg / delta_time_C

      ! Reduce the above ground yasso pools accordingly
      fraction = 1._lp - burned_fract(ic)
      c_acid_ag1      (ic) = fraction * c_acid_ag1       (ic)
      c_water_ag1     (ic) = fraction * c_water_ag1      (ic)
      c_ethanol_ag1   (ic) = fraction * c_ethanol_ag1    (ic)
      c_nonsoluble_ag1(ic) = fraction * c_nonsoluble_ag1 (ic)

      c_acid_ag2      (ic) = fraction * c_acid_ag2       (ic)
      c_water_ag2     (ic) = fraction * c_water_ag2      (ic)
      c_ethanol_ag2   (ic) = fraction * c_ethanol_ag2    (ic)
      c_nonsoluble_ag2(ic) = fraction * c_nonsoluble_ag2 (ic)

      ! Add below ground fraction of the burned vegetation to the below ground green litter pools
      fraction = (1._lp - fract_green_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the remaings of burned wood to the yasso litter pools
      !  Above ground
      fraction = (1._lp - fire_fract_wood_2_atmos) * fract_wood_aboveGround * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._lp - fract_wood_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._lp - burned_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_fire


  !> This routine handles the impact of fire on the carbon pools.
  !>
  !> We assume that the above ground litter carbon (yasso ag1 and ag2 pools), is completely released
  !> to the atmosphere while the carbon of the living plants (c_green, c_reserve, c_woods) is partly
  !> released to the atmosphere and partly relocated into the litter pools.
  !
  SUBROUTINE relocate_carbon_vertical_fire( &
    & delta_time_C,                      & ! in
    & burned_fract,                      & ! in
    & fire_fract_wood_2_atmos,           & ! in
    & LeafLit_coef,                      & ! in
    & WoodLit_coef,                      & ! in
    & fuel_moisture,                     & ! in
    & c_green,                           & ! inout
    & c_woods,                           & ! inout
    & c_reserve,                         & ! inout
    & cflux_dist_green_2_soil,           & ! inout
    & cflux_dist_woods_2_soil,           & ! inout
    & c_acid_ag1,                        & ! inout
    & c_water_ag1,                       & ! inout
    & c_ethanol_ag1,                     & ! inout
    & c_nonsoluble_ag1,                  & ! inout
    & c_acid_ag2,                        & ! inout
    & c_water_ag2,                       & ! inout
    & c_ethanol_ag2,                     & ! inout
    & c_nonsoluble_ag2,                  & ! inout
    & c_acid_bg1,                        & ! inout
    & c_water_bg1,                       & ! inout
    & c_ethanol_bg1,                     & ! inout
    & c_nonsoluble_bg1,                  & ! inout
    & c_acid_bg2,                        & ! inout
    & c_water_bg2,                       & ! inout
    & c_ethanol_bg2,                     & ! inout
    & c_nonsoluble_bg2,                  & ! inout
    & c_humus_1,                         & ! inout
    & c_humus_2,                         & ! inout
    & c_acid_bg1_sl,                     & ! inout
    & c_water_bg1_sl,                    & ! inout
    & c_ethanol_bg1_sl,                  & ! inout
    & c_nonsoluble_bg1_sl,               & ! inout
    & c_acid_bg2_sl,                     & ! inout
    & c_water_bg2_sl,                    & ! inout
    & c_ethanol_bg2_sl,                  & ! inout
    & c_nonsoluble_bg2_sl,               & ! inout
    & c_humus_1_sl,                      & ! inout
    & c_humus_2_sl,                      & ! inout
    & co2flux_fire_all_2_atm             & ! out
    & )

    USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround, molarMassCO2_kg

    ! Arguments
    REAL(lp), INTENT(in)    :: delta_time_C                 !< Carbon time step [s]
    REAL(lp), INTENT(in)    :: burned_fract(:)              !< Fraction of the vegetated area of each tile burned till the last
                                                            !< call of this routine
    REAL(lp), INTENT(in)    :: fire_fract_wood_2_atmos      !< Fraction of above ground wood immediately emitted to atmosphere
                                                            !< in case of fire
    REAL(lp), INTENT(in)    :: fuel_moisture(:)             !
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)              !< Factor to spread non woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)              !< Factor to spread woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(inout) :: c_green(:)                   !< Amount of green carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_woods(:)                   !< Amount of woody carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_reserve(:)                 !< Amount of reserve carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: cflux_dist_green_2_soil(:)   !< Carbon relocated by disturbances from non-woody biomass
                                                            !< to the green litter pools [mol(C)/m^2(canopy)s-1]
    REAL(lp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   !< Carbon relocated by disturbances from woody biomass to
                                                            !< the wood litter pools [mol(C)/m^2(canopy)s-1]
    !   YASSO
    ! SIZE CLASS 1: Green, originating from non-woody biomass
    !   above ground C pools
    REAL(lp), INTENT(inout) :: c_acid_ag1(:)          !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_ag1(:)         !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_ag1(:)       !< Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1(:)    !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp), INTENT(inout) :: c_acid_bg1(:)          !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg1(:)         !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1(:)       !< Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1(:)    !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_1(:)           !< Humus pool from green litter [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2: Woody, originating from woody biomass
    !   above ground C pools
    REAL(lp), INTENT(inout) :: c_acid_ag2(:)          !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_ag2(:)         !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_ag2(:)       !< Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2(:)    !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp), INTENT(inout) :: c_acid_bg2(:)          !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg2(:)         !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2(:)       !< Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2(:)    !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_2(:)           !< Humus pool from woody litter [mol(C)/m^2(canopy)]

    REAL(lp), INTENT(inout) :: c_acid_bg1_sl(:,:)          !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg1_sl(:,:)         !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1_sl(:,:)       !< Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1_sl(:,:)    !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_1_sl(:,:)           !< Humus pool from green litter [mol(C)/m^2(canopy)]

    REAL(lp), INTENT(inout) :: c_acid_bg2_sl(:,:)          !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg2_sl(:,:)         !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2_sl(:,:)       !< Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2_sl(:,:)    !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_2_sl(:,:)           !< Humus pool from woody litter [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(out)   :: co2flux_fire_all_2_atm(:)   !< CO2 released by fire [kg(CO2)m-2(canopy)s-1]

    ! Local variables
    REAL(lp)  :: non_woody_litter, fraction ! Helper
    INTEGER   :: ic, nc

    nc = SIZE(burned_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc
      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Carbon released due to disturbances from the green and reserve pools, that enters the green litter and
      ! soil pools. We here add the contribution of fires.
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + &
        & (c_green(ic) + c_reserve(ic)) * (1._lp - fract_green_aboveGround) * burned_fract(ic) / delta_time_C
      ! Carbon released due to disturbances from the woody vegetation pool, that enters the woody litter and soil
      ! pools. We here add the contribution of fires.
      cflux_dist_woods_2_soil(ic)  = cflux_dist_woods_2_soil(ic) + &
        & c_woods(ic) * (1._lp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * burned_fract(ic) / delta_time_C

      ! Determine the carbon flux to the atmosphere due to fire: We assume that the above ground non-woody
      ! biomass and above ground litter are burnt completely, while parts of the above ground woody biomass
      ! and below ground matter resists the flames.
      co2flux_fire_all_2_atm(ic) = ( &
          ! above ground fraction of living biomass
        & non_woody_litter * fract_green_aboveGround                          &
          ! add the above ground YASSO pools
        & + c_water_ag1(ic) + c_acid_ag1(ic) + c_ethanol_ag1(ic)              &
        & + c_nonsoluble_ag1(ic)                                              &
        & + c_water_ag2(ic) + c_acid_ag2(ic) + c_ethanol_ag2(ic)              &
        & + c_nonsoluble_ag2(ic)                                              &
          ! add the fire affected fraction of the above ground wood pools
        & + c_woods(ic) * fire_fract_wood_2_atmos * fract_wood_aboveGround    &
          ! scale with the burned fraction
        & ) * burned_fract(ic) * molarMassCO2_kg / delta_time_C

      ! Reduce the above ground yasso pools accordingly
      fraction = 1._lp - burned_fract(ic)
      c_acid_ag1      (ic) = fraction * c_acid_ag1       (ic)
      c_water_ag1     (ic) = fraction * c_water_ag1      (ic)
      c_ethanol_ag1   (ic) = fraction * c_ethanol_ag1    (ic)
      c_nonsoluble_ag1(ic) = fraction * c_nonsoluble_ag1 (ic)

      c_acid_ag2      (ic) = fraction * c_acid_ag2       (ic)
      c_water_ag2     (ic) = fraction * c_water_ag2      (ic)
      c_ethanol_ag2   (ic) = fraction * c_ethanol_ag2    (ic)
      c_nonsoluble_ag2(ic) = fraction * c_nonsoluble_ag2 (ic)

      ! Add the below ground fraction of the burned vegetation to the below ground green litter pools
      fraction = (1._lp-fract_green_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1_sl(ic,1), c_water_bg1_sl(ic,1), c_ethanol_bg1_sl(ic,1), &
        &                          c_nonsoluble_bg1_sl(ic,1), c_humus_1_sl(ic,1), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the remanings of burned wood to the yasso litter pools
      !  Above ground
      fraction = (1._lp-fire_fract_wood_2_atmos) * fract_wood_aboveGround * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._lp - fract_wood_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2_sl(ic,1), c_water_bg2_sl(ic,1), c_ethanol_bg2_sl(ic,1), &
        &                          c_nonsoluble_bg2_sl(ic,1), c_humus_2_sl(ic,1), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._lp - burned_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction

      c_acid_bg1(ic)       = SUM(c_acid_bg1_sl(ic,:))
      c_water_bg1(ic)      = SUM(c_water_bg1_sl(ic,:))
      c_ethanol_bg1(ic)    = SUM(c_ethanol_bg1_sl(ic,:))
      c_nonsoluble_bg1(ic) = SUM(c_nonsoluble_bg1_sl(ic,:))
      c_humus_1(ic)        = SUM(c_humus_1_sl(ic,:))
      c_acid_bg2(ic)       = SUM(c_acid_bg2_sl(ic,:))
      c_water_bg2(ic)      = SUM(c_water_bg2_sl(ic,:))
      c_ethanol_bg2(ic)    = SUM(c_ethanol_bg2_sl(ic,:))
      c_nonsoluble_bg2(ic) = SUM(c_nonsoluble_bg2_sl(ic,:))
      c_humus_2(ic)        = SUM(c_humus_2_sl(ic,:))
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_vertical_fire


  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (c_green_st, c_reserve_st and
  ! c_woods_st) is partly relocated into the litter pools.
  !
  SUBROUTINE relocate_carbon_damage(  &
    & delta_time_C,                   & ! in
    & damaged_fract,                  & ! in
    & LeafLit_coef,                   & ! in
    & WoodLit_coef,                   & ! in
    & c_green,                        & ! inout
    & c_woods,                        & ! inout
    & c_reserve,                      & ! inout
    & cflux_dist_green_2_soil,        & ! inout
    & cflux_dist_woods_2_soil,        & ! inout
    & c_acid_ag1,                     & ! inout
    & c_water_ag1,                    & ! inout
    & c_ethanol_ag1,                  & ! inout
    & c_nonsoluble_ag1,               & ! inout
    & c_acid_ag2,                     & ! inout
    & c_water_ag2,                    & ! inout
    & c_ethanol_ag2,                  & ! inout
    & c_nonsoluble_ag2,               & ! inout
    & c_acid_bg1,                     & ! inout
    & c_water_bg1,                    & ! inout
    & c_ethanol_bg1,                  & ! inout
    & c_nonsoluble_bg1,               & ! inout
    & c_acid_bg2,                     & ! inout
    & c_water_bg2,                    & ! inout
    & c_ethanol_bg2,                  & ! inout
    & c_nonsoluble_bg2,               & ! inout
    & c_humus_1,                      & ! inout
    & c_humus_2                       & ! inout
    & )

    ! Arguments
    REAL(lp), INTENT(in)   :: delta_time_C             !< Carbon time step [s]
    REAL(lp),INTENT(in)    :: damaged_fract(:)         ! Fraction of the vegetated area of each tile damaged since the last
                                                       ! call of this routine
    REAL(lp),INTENT(in)    :: LeafLit_coef(:)          ! Factor to spread non woody litter into yasso pools [ ]
    REAL(lp),INTENT(in)    :: WoodLit_coef(:)          ! Factor to spread woody litterfall into yasso pools [ ]
    REAL(lp),INTENT(inout) :: c_green(:)               ! Value of green carbon pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_woods(:)               ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_reserve(:)             ! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    REAL(lp),INTENT(inout) :: cflux_dist_green_2_soil(:)          ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
    REAL(lp),INTENT(inout) :: cflux_dist_woods_2_soil(:)          ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
                                                                  ! .. to the green litter pools [mol(C)/m^2(canopy)s-1]

    !   YASSO
    ! Size class 1: green
    !   above ground C pools
    REAL(lp),INTENT(inout) :: c_acid_ag1(:)        ! Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_ag1(:)       ! Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_ag1(:)     ! Ethanol soluble green litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_ag1(:)  ! Non-soluble green litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp),INTENT(inout) :: c_acid_bg1(:)        ! Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_bg1(:)       ! Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_bg1(:)     ! Ethanol soluble green litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_bg1(:)  ! Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_humus_1(:)         ! Humus pool from green litter [mol(C)/m^2(canopy)]
    ! Size class 2: wood
    !   above ground C pools
    REAL(lp),INTENT(inout) :: c_acid_ag2(:)        ! Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_ag2(:)       ! Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_ag2(:)     ! Ethanol soluble woody litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_ag2(:)  ! Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp),INTENT(inout) :: c_acid_bg2(:)        ! Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_bg2(:)       ! Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_bg2(:)     ! Ethanol soluble woody litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_bg2(:)  ! Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_humus_2(:)         ! Humus pool from woody litter [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(lp)  :: non_woody_litter, fraction  ! Helper
    INTEGER   :: ic, nc

    ! Initializations
    nc = SIZE(damaged_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc

      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Amount of carbon transfered from damaged wood pool to the wood litter pools
      cflux_dist_woods_2_soil(ic) = cflux_dist_woods_2_soil(ic) + c_woods(ic) * damaged_fract(ic) / delta_time_C

      ! Amount of carbon transfered from the damaged non-woody vegetation pools to the green litter pools
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + non_woody_litter * damaged_fract(ic) / delta_time_C

      ! Add the carbon of damaged non-woody vegetation pools to the yasso soil pools
      !  Above ground
      fraction = fract_green_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), c_nonsoluble_ag1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)
      !  Below ground
      fraction = (1._lp - fract_green_aboveGround) * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Add the carbon of damaged woody vegetation pool to the yasso soil pools
      !  Above ground
      fraction = fract_wood_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._lp - fract_wood_aboveGround)  * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._lp - damaged_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_damage

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  !< In this routine we model the impact of vegetation damage (i.e. wind through) on the carbon pools:
  !< Depending on the disturbed area fraction carbon from the living vegetation pools (c_green,
  !< c_reserve and c_woods) is relocated  to the litter and soil pools.
  !
  SUBROUTINE relocate_carbon_vertical_damage(  &
    & delta_time_C,                   & ! in
    & damaged_fract,                  & ! in
    & LeafLit_coef,                   & ! in
    & WoodLit_coef,                   & ! in
    & c_green,                        & ! inout
    & c_woods,                        & ! inout
    & c_reserve,                      & ! inout
    & cflux_dist_green_2_soil,        & ! inout
    & cflux_dist_woods_2_soil,        & ! inout
    & c_acid_ag1,                     & ! inout
    & c_water_ag1,                    & ! inout
    & c_ethanol_ag1,                  & ! inout
    & c_nonsoluble_ag1,               & ! inout
    & c_acid_ag2,                     & ! inout
    & c_water_ag2,                    & ! inout
    & c_ethanol_ag2,                  & ! inout
    & c_nonsoluble_ag2,               & ! inout
    & c_acid_bg1,                     & ! inout
    & c_water_bg1,                    & ! inout
    & c_ethanol_bg1,                  & ! inout
    & c_nonsoluble_bg1,               & ! inout
    & c_acid_bg2,                     & ! inout
    & c_water_bg2,                    & ! inout
    & c_ethanol_bg2,                  & ! inout
    & c_nonsoluble_bg2,               & ! inout
    & c_humus_1,                      & ! inout
    & c_humus_2,                      & ! inout
    & c_acid_bg1_sl,                  & ! inout
    & c_water_bg1_sl,                 & ! inout
    & c_ethanol_bg1_sl,               & ! inout
    & c_nonsoluble_bg1_sl,            & ! inout
    & c_acid_bg2_sl,                  & ! inout
    & c_water_bg2_sl,                 & ! inout
    & c_ethanol_bg2_sl,               & ! inout
    & c_nonsoluble_bg2_sl,            & ! inout
    & c_humus_1_sl,                   & ! inout
    & c_humus_2_sl                    & ! inout
    & )

    USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround

    ! Arguments
    REAL(lp), INTENT(in)    :: delta_time_C                !< Carbon time step [s]
    REAL(lp), INTENT(in)    :: damaged_fract(:)            !< Fraction of the vegetated area of each tile damaged since the last
                                                           !< call of this routine
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)             !< Factor to spread non woody litter into yasso pools [ ]
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)             !< Factor to spread woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(inout) :: c_green(:)                  !< Amount of carbon in green pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_woods(:)                  !< Amount carbon in wood pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_reserve(:)                !< Amount of carbon in reserve pool [mol(C)/m^(canopy)2]
    REAL(lp), INTENT(inout) :: cflux_dist_green_2_soil(:)  !< Disturbance-induced carbon flux of non-woody living biomass
                                                           !< to the green litter and soil pools [mol(C)/m^2(canopy)s-1]
    REAL(lp), INTENT(inout) :: cflux_dist_woods_2_soil(:)  !< Disturbance-induced carbon flux of woody living biomass to
                                                           !< the woody litter and soil pools [mol(C)/m^2(canopy)s-1]

    !   YASSO
    ! Size class 1: green
    !   above ground C pools
    REAL(lp),INTENT(inout) :: c_acid_ag1(:)        !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_ag1(:)       !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_ag1(:)     !< Ethanol soluble green litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_ag1(:)  !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp),INTENT(inout) :: c_acid_bg1(:)        !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_bg1(:)       !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_bg1(:)     !< Ethanol soluble green litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_bg1(:)  !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_humus_1(:)         !< Humus pool from green litter [mol(C)/m^2(canopy)]
    ! Size class 2: wood
    !   above ground C pools
    REAL(lp),INTENT(inout) :: c_acid_ag2(:)        !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_ag2(:)       !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_ag2(:)     !< Ethanol soluble woody litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_ag2(:)  !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp),INTENT(inout) :: c_acid_bg2(:)        !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_bg2(:)       !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_bg2(:)     !< Ethanol soluble woody litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_bg2(:)  !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_humus_2(:)         !< Humus pool from woody litter [mol(C)/m^2(canopy)]

    REAL(lp), INTENT(inout) :: c_acid_bg1_sl(:,:)          !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg1_sl(:,:)         !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1_sl(:,:)       !< Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1_sl(:,:)    !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_1_sl(:,:)           !< Humus pool from green litter [mol(C)/m^2(canopy)]

    REAL(lp), INTENT(inout) :: c_acid_bg2_sl(:,:)          !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg2_sl(:,:)         !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2_sl(:,:)       !< Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2_sl(:,:)    !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_2_sl(:,:)           !< Humus pool from woody litter [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(lp)  :: non_woody_litter, fraction  ! Helper
    INTEGER   :: ic, nc

    ! Initializations
    nc = SIZE(damaged_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc

      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Carbon transferred from the wood pool to the woody litter and soil pools due to disturbances;
      ! we here add the contribution of wind throw.
      cflux_dist_woods_2_soil(ic) = cflux_dist_woods_2_soil(ic) + c_woods(ic) * damaged_fract(ic) / delta_time_C

      ! Carbon transferred from the non-woody vegetation pools to the green litter and soil pools due to
      ! disturbances; we here add the contribution of wind throw.
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + non_woody_litter * damaged_fract(ic) / delta_time_C

      ! Distribute the carbon of damaged non-woody vegetation to the yasso litter and soil pools
      !  Above ground
      fraction = fract_green_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), c_nonsoluble_ag1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)
      !  Below ground
      fraction = (1._lp - fract_green_aboveGround) * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1_sl(ic,1), c_water_bg1_sl(ic,1), c_ethanol_bg1_sl(ic,1), &
        &                          c_nonsoluble_bg1_sl(ic,1), c_humus_1_sl(ic,1), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the carbon of damaged woody vegetation to the yasso litter and soil pools
      !  Above ground
      fraction = fract_wood_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._lp - fract_wood_aboveGround)  * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2_sl(ic,1), c_water_bg2_sl(ic,1), c_ethanol_bg2_sl(ic,1), &
        &                          c_nonsoluble_bg2_sl(ic,1), c_humus_2_sl(ic,1), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._lp - damaged_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction

      c_acid_bg1(ic)       = SUM(c_acid_bg1_sl(ic,:))
      c_water_bg1(ic)      = SUM(c_water_bg1_sl(ic,:))
      c_ethanol_bg1(ic)    = SUM(c_ethanol_bg1_sl(ic,:))
      c_nonsoluble_bg1(ic) = SUM(c_nonsoluble_bg1_sl(ic,:))
      c_humus_1(ic)        = SUM(c_humus_1_sl(ic,:))
      c_acid_bg2(ic)       = SUM(c_acid_bg2_sl(ic,:))
      c_water_bg2(ic)      = SUM(c_water_bg2_sl(ic,:))
      c_ethanol_bg2(ic)    = SUM(c_ethanol_bg2_sl(ic,:))
      c_nonsoluble_bg2(ic) = SUM(c_nonsoluble_bg2_sl(ic,:))
      c_humus_2(ic)        = SUM(c_humus_2_sl(ic,:))

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_vertical_damage


#endif
END MODULE mo_carbon_process
