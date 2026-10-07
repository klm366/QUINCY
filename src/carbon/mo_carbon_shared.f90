!> Contains carbon process routines shared between the standard carbon usage and the methane version
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
MODULE mo_carbon_shared
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_plant_c, yasso, distribute_yasso_litter, add_litter_to_yasso_pool, get_per_tile

  CHARACTER(*), PARAMETER :: modname = 'mo_carbon_shared'

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
    REAL(wp), INTENT(in)    :: delta_time_C            !< Carbon time step [s]
    LOGICAL,  INTENT(in)    :: is_crop                 !< Is this a crop tile?
    REAL(wp), INTENT(in)    :: LAI                     !< Yesterdays mean LAI
    REAL(wp), INTENT(in)    :: LAI_previous            !< The day before yesterdays mean LAI
    REAL(wp), INTENT(in)    :: max_LAI                 !< Maximum value of LAI
    REAL(wp), INTENT(in)    :: NPP_pot_yDayMean        !< Yesterdays mean NPP rate
    REAL(wp), INTENT(in)    :: fract_npp_2_woodPool    !< Fraction of NPP to be put maximally into the green pool
    REAL(wp), INTENT(in)    :: fract_npp_2_reservePool !< Fraction of NPP to be put into the optimally into the reserve pool
    REAL(wp), INTENT(in)    :: fract_npp_2_exudates    !< Fraction of NPP to be put into the optimally into the root exudates
    REAL(wp), INTENT(in)    :: fract_green_2_herbivory
    REAL(wp), INTENT(in)    :: tau_c_woods             !< Time constant by which woods Pool is depreciated  [d]
                                                       !< Note: in lctlib it is given in [yr-1] but used by JS4 in [d-1].
    REAL(wp), INTENT(in)    :: LAI_shed_constant       !< Leaf shedding at a constant rate for evergreens etc. [d-1]
    REAL(wp), INTENT(in)    :: Max_C_content_woods     !< Maximum carbon content of wood pool (from lctLib)
    REAL(wp), INTENT(in)    :: specific_leaf_area_C    !< Specific leaf area (from lctLib) [m2(leaf) mol(C)-1]
    REAL(wp), INTENT(in)    :: reserveC2leafC          !< Ratio of max. carbon in reserve pool to max. carbon of leaves
    ! Plant carbon pools: on input: last value; updated on output
    REAL(wp), INTENT(inout) :: c_green                 !< Green carbon pool
    REAL(wp), INTENT(inout) :: c_woods                 !< Wood carbon pool
    REAL(wp), INTENT(inout) :: c_reserve               !< Reserve carbon pool
    REAL(wp), INTENT(inout) :: c_crop_harvest          !< Crop harvest carbon pool
    ! Output arguments
    REAL(wp), INTENT(out)   :: NPP_flux_correction     !< Amount by which the NPP rate entering the routine has to be corrected. This
                                                       !<   correction arises either because otherwise the reserve pool would get
                                                       !<   negative (positive correction), or the wood pool would exceed its
                                                       !<   maximum value (negative correction).
    REAL(wp), INTENT(out)   :: NPP_excess              !< Part of NPP that could not be stored in one of the plant carbon pools
                                                       !<   (green, wood, reserve), but had to be  thrown away (a posteriori
                                                       !<   reduction of NPP)
    REAL(wp), INTENT(out)   :: NPP_act                 !< Actual NPP after N-limitation and excess carbon drop
    REAL(wp), INTENT(out)   :: Cflx_2_root_exudates    !< Carbon flux into root exudates
    REAL(wp), INTENT(out)   :: Cflx_greenwood_2_litter !< Total carbon flux from vegetation to litter (green+woody, ag+bg)
    REAL(wp), INTENT(out)   :: Cflx_green_2_litter     !< Total carbon flux of easily decomposable material
    REAL(wp), INTENT(out)   :: Cflx_wood_2_litter      !< Total carbon flux of slowly decomposable material
    REAL(wp), INTENT(out)   :: Cflx_herbivory
    REAL(wp), INTENT(out)   :: Cflx_herb_2_littergreen
    REAL(wp), INTENT(out)   :: Cflx_herb_2_atm
    REAL(wp), INTENT(out)   :: Cflx_crop_harvest_2_atm

    ! Parameters
    REAL(wp), PARAMETER :: days_per_year      = 365.0_wp   !< Just for testing - should be imported, but slightly different value
    ! R: Normally parameters should be in mo_carbon_constants, but the following parameters are used only in this subroutine!
    !    Therefore I put them here. Is there any not structural advantage to put them into mo_carbon_constants?
    REAL(wp), PARAMETER :: greenC2leafC         = 4.0_wp          !< ratio of carbon green pool (leaves, fine roots, starches,
                                                                  !<   sugars) to leaf carbon
    REAL(wp), PARAMETER :: fract_c_crop_harvest = 0.5_wp          !< fraction of C from ag litter flux put to the crop harvest pool
    REAL(wp), PARAMETER :: fract_c_faeces2_LG   = 0.3_wp          !< fraction of C from faces put into the litter green pool
    REAL(wp), PARAMETER :: tau_c_reserve        = 1.0_wp * days_per_year !< time constant by which reserve pool is depreciated [yr]
    REAL(wp), PARAMETER :: tau_c_crop_harvest   = 1.0_wp * days_per_year !< time constant by which the crop harvest pool decays [yr]

    REAL(wp) :: c_green_max, c_reserve_optimal
    REAL(wp) :: c_NPP_2_greenPool, c_NPP_2_woodPool, c_NPP_2_reservePool, c_NPP_2_rootExudates, c_green_2_herbivory
    REAL(wp) :: c_2_litter_greenPool, c_2_litter_woodPool, c_reserve_2_litter_greenPool
    REAL(wp) :: c_excess, c_crop_harvest_2_atm
    REAL(wp) :: leaf_shedding

    REAL(wp) :: c_green_pot, c_woods_pot         !! values of potential carbon pools i.e. before accounting for N-limitation
    REAL(wp) :: c_reserve_pot

    REAL(wp) :: c_NPP_2_green_pot, c_NPP_2_wood_pot, c_NPP_2_reserve_pot
    REAL(wp) :: c_2_crop_harvest
    REAL(wp) :: c_faeces_2_LG, c_faeces_2_atm
    ! Preparation for N-limitation
    REAL(wp) :: c_NPP_2_green_via_mobileN, c_NPP_2_wood_via_mobileN
    REAL(wp) :: redFact_Nlimitation
    ! R: only for spitfire
    !REAL(wp) :: litter_wood_new, litter_wood_total

    ! Initializations
    c_2_crop_harvest     = 0.0_wp

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
       MIN(greenC2leafC * leaf_shedding / specific_leaf_area_c, &   !!    by leaf and root shedding to the green litter pool
           c_green)
    c_green = c_green - c_2_litter_greenPool                        !! This removes carbon from the green pool
    c_excess = MAX(0.0_wp,c_green - c_green_max)                    !! If by the reduction of LAI the maximum value of the green
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
    c_faeces_2_atm      = c_green_2_herbivory * (1.0_wp - fract_c_faeces2_LG) !! c_faeces_2_atm to be lost to atmosphere
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

    NPP_flux_correction = 0.0_wp
    IF (NPP_pot_yDayMean <= 0.0_wp) THEN                      !! In case of negative or zero NPP
      NPP_excess = 0.0_wp

      c_NPP_2_reservePool =  NPP_pot_yDayMean * delta_time_C  !! .. the plant looses C and
      c_reserve     = c_reserve + c_NPP_2_reservePool         !! .. we try to take it from the reserve pool
      IF (c_reserve < 0.0_wp) THEN                            !! But if thereby the reserve pool gets negative
        c_NPP_2_reservePool = c_NPP_2_reservePool - c_reserve !! .. we absorb negative NPP only up to zero reserve pool and
        NPP_flux_correction = -c_reserve / delta_time_C       !! .. give the deficit as a flux correction back to the
        c_reserve = 0.0_wp                                    !! .. calling routine. Hence all Carbon from the reserve pool
      END IF                                                  !! .. is used up.

      !! For correct handling of nitrogen limitation below (which handles only the case of positive NPP) set NPP to zero:
      c_NPP_2_reservePool  = 0.0_wp
      c_NPP_2_greenPool    = 0.0_wp    !! All other transfer rates are zero
      c_NPP_2_woodPool     = 0.0_wp
      Cflx_2_root_exudates = 0.0_wp
      c_reserve_pot        = c_reserve !! remember modified reserve pool

    ELSE  !! NPP is positive

      !! NPP it is distributed according to predefined relative fraction
      c_NPP_2_woodPool     = delta_time_C * fract_npp_2_woodPool    * NPP_pot_yDayMean
      c_NPP_2_reservePool  = delta_time_C * fract_npp_2_reservePool * NPP_pot_yDayMean
      c_NPP_2_rootExudates = delta_time_C * fract_npp_2_exudates    * NPP_pot_yDayMean
      c_NPP_2_greenPool    = delta_time_C * &
        & (1.0_wp-fract_npp_2_woodPool-fract_npp_2_reservePool-fract_npp_2_exudates) * NPP_pot_yDayMean

      !! Growth of Wood Pool (structural carbon of living plants)

      c_woods_pot       = c_woods_pot + c_NPP_2_woodPool                !! Then it is attempted to put all NPP share into it
      c_excess          = MAX(0.0_wp,c_woods_pot-Max_C_content_woods)   !! .. but thereby the pool may get too large
      c_woods_pot       = c_woods_pot - c_excess                        !! .. so that the excess carbon has once more to be
                                                                        !! .. subtracted
      c_NPP_2_greenPool = c_NPP_2_greenPool + c_excess                  !! .. and instead made available to the green pool
      c_NPP_2_woodPool  = c_NPP_2_woodPool  - c_excess                  !! .. and the actual amount of NPP put to the wood pool
                                                                        !! .. is less

      !! Growth of Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
      IF (c_reserve_pot < c_reserve_optimal) THEN                       !! .. If the reserve pool is smaller than optimal,
        c_reserve_pot     = c_reserve_pot     + c_NPP_2_reservePool     !! .... it is filled by the available NPP
        c_excess          =                                         &   !! .... Thereby it may happen that it gets larger than
             MAX(0.0_wp,c_reserve_pot - c_reserve_optimal)              !! .... optimal so that there is excess carbon
        c_reserve_pot       = c_reserve_pot       - c_excess            !! .... that needs not be taken up
        c_NPP_2_greenPool   = c_NPP_2_greenPool   + c_excess            !! .... but can better be used to increase the green
        c_NPP_2_reservePool = c_NPP_2_reservePool - c_excess            !! .... pool and the actual amount of NPP put to the
                                                                        !! .... reserve pool is less.
      ELSE                                                              !! .... Otherwise (reserve pool is larger than optimal)
        c_NPP_2_greenPool   = c_NPP_2_greenPool   + c_NPP_2_reservePool !! .... all NPP is left for the green pool
        c_NPP_2_reservePool = 0.0_wp                                    !! .... so that nothing is stored in the reserve pool.
      END IF

      !! Growth of Green Pool (leaves and fine roots): (in case of too much NPP, try to put it into the reserve pool).

      c_green_pot       = c_green_pot + c_NPP_2_greenPool               !! Green pool is filled by the available NPP.
      c_excess          = MAX(0.0_wp,c_green_pot - c_green_max)         !! .. Thereby it may get larger than appropriate for
                                                                        !!    current LAI.
      c_NPP_2_greenPool = c_NPP_2_greenPool - c_excess                  !! .. Hence the actual amount of NPP put to the green
                                                                        !!    pool is less.
      c_green_pot       = c_green_pot - c_excess                        !! .. This excess carbon needs not be taken up but
      IF (c_reserve_pot < c_reserve_optimal) THEN                       !! .... if the reserve pool is smaller than optimal
        c_reserve_pot       = c_reserve_pot       + c_excess            !! .... it is tried to put the carbon there,
        c_NPP_2_reservePool = c_NPP_2_reservePool + c_excess            !! .... which means that additional NPP is put to the
                                                                        !!      reserve pool.
        c_excess            =                                   &       !! .... Thereby it may happen that the reserve pool
                 MAX(0.0_wp,c_reserve_pot         - c_reserve_optimal)  !!      increases beyond the optimal value.
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

    c_NPP_2_green_via_mobileN = 0.0_wp
    c_NPP_2_wood_via_mobileN  = 0.0_wp

    c_NPP_2_green_pot = c_NPP_2_green_pot - c_NPP_2_green_via_mobileN
    c_NPP_2_wood_pot  = c_NPP_2_wood_pot  - c_NPP_2_wood_via_mobileN

    ! R: These lines were calculated also with_Nitrogen=.FALSE. in JS3. However, it was not used further when with_Nitrogen=.FALSE.
    !    In JS3 it should stand within the brackets.
    !Nflux_NPP_2_green_pot = c_NPP_2_green_pot/cn_green
    !Nflux_NPP_2_wood_pot  = c_NPP_2_wood_pot/cn_woods
    !minN_plant_demand = Nflux_NPP_2_green_pot + Nflux_NPP_2_wood_pot

    !! ============== COLLECT POTENTIAL C-FLUXES ===============================================================

    c_NPP_2_reserve_pot = c_NPP_2_reservePool !! potential growth of reserve carbon
    redFact_Nlimitation = 1.0_wp              !! This assures absence of nitrogen limitation in this mixed C-N code

    !! ============== UPDATE ALL C-POOLS =======================================================================

    c_green   = c_green   + c_NPP_2_green_pot * redFact_Nlimitation
    c_woods   = c_woods   + c_NPP_2_wood_pot  * redFact_Nlimitation
    c_reserve = c_reserve + c_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N

    !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

    IF (NPP_pot_yDayMean > 0.0_wp) THEN
      NPP_act = (c_NPP_2_green_via_mobileN + c_NPP_2_wood_via_mobileN + c_NPP_2_rootExudates +             &
  & redFact_Nlimitation * (c_NPP_2_green_pot + c_NPP_2_wood_pot) + c_NPP_2_reserve_pot  &
  & ) / delta_time_C
    ELSE
      NPP_act = NPP_pot_yDayMean + NPP_flux_correction  ! NPP_pot_yDayMean is negative: NPP_flux_correction is zero
                                                         ! except c_reserve was < 0.0_wp.
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
    !           if ((litter_wood_total >  EPSILON(1._wp)) .AND. (litter_wood_new < litter_wood_total)) then
    !              fract_litter_wood_new = litter_wood_new / litter_wood_total
    !           else
    !              fract_litter_wood_new = 1._wp    ! fract_litter_wood_new=1 means fuel fractions will be reset to original values
    !           endif
    !        endif


  END SUBROUTINE calc_plant_c

#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE yasso (delta_time_C, Yasso_in_pools, Weather, c_2_litter, Lit_coefV, WoodLitterSize, Yasso_out, fract_aboveground, &
  &               c_2_rootExudates, redFact_Nlimit, is_anoxic, is_peat,                                             &
  &               A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H)
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

    USE mo_carbon_constants, ONLY: sec_per_year, LCT_HUMUS,                &
  & POOL_ACID_AG, POOL_WATER_AG, POOL_ETHANOL_AG,         &
  & POOL_NONSOLUBLE_AG,                                   &
  & POOL_ACID_BG, POOL_WATER_BG, POOL_ETHANOL_BG,         &
  & POOL_NONSOLUBLE_BG, POOL_HUMUS,                       &
  & FLUX_RESP_TOTAL, FLUX_2_HUMUS_AG,                     &
  & FLUX_2_HUMUS_BG, FLUX_FROM_HUMUS_DIAG,                &
  & FLUX_RESP_LITTER_AG, FLUX_RESP_LITTER_BG,             &
  & REL_FLUX_AWEN, FLUX_FROM_HUMUS, FLUX_LITTER_2_HUMUS
    USE mo_jsb_physical_constants, ONLY: tmelt
    IMPLICIT NONE

    REAL(wp),                INTENT(in)  :: delta_time_C       !< Carbon time step [s]
    ! TODO: instead of Yasso_in_pools and Yasso_out use a variable Yasso_pools or such with intent INOUT
    !       for the Yasso pools and a variable Yasso_fluxes or such replacing Yasso_out(10:18) with intent out.

    REAL(wp), DIMENSION(9),  INTENT(in)  :: Yasso_in_pools     !< Yasso pools in [mol(C) m-2]
    REAL(wp),                INTENT(in)  :: WoodLitterSize     !< Size of woody litter;
                                                               !<   is zero when the subroutine is called for non-woody litter
    REAL(wp), DIMENSION(2),  INTENT(in)  :: Weather            !< climatic drivers: air temperature and precipitation, 15-d means
                                                               !< Weather(1) in [K] ; Weather(2) in [mm s-1]
    REAL(wp),                INTENT(in)  :: c_2_litter         !< fresh litter inputs  (above and belowground) [mol(C) m-2 d-1]
    REAL(wp), DIMENSION(5),  INTENT(in)  :: Lit_coefV          !< fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(wp),                INTENT(in)  :: fract_aboveground  !< parameter to seperate litter influx into above- and belowground
                                                               !<     part                                                   [ ]
    REAL(wp),                INTENT(in)  :: c_2_rootExudates   !< root exudates flux [mol(C) m-2 d-1]

    REAL(wp), DIMENSION(18), INTENT(out) :: Yasso_out          !< updated Y pools (1:9) and fluxes (10:18)
                                                               !<   see update_Cpools for more information

    REAL(wp), OPTIONAL,      INTENT(in)  :: redFact_Nlimit
    LOGICAL,  OPTIONAL,      INTENT(in)  :: is_anoxic
    LOGICAL,  OPTIONAL,      INTENT(in)  :: is_peat

    ! Possibilities to overwrite (some of) the default exchange coefficients (AWEN2H for backward compatibility)
    REAL(wp), OPTIONAL,      INTENT(in)  :: A2W, A2E, A2N, W2A, W2E, W2N, E2A, E2W, E2N, N2A, N2W, N2E, AWEN2H, N2H, AWE2H

    ! Yasso parameters (Tuomi '11)
    !   reference decompositation rates for T = 0C and unlimited water availability in [yr-1]...
    REAL(wp), PARAMETER   :: ref_decomp_rate_acid       = -0.72_wp    !<  ... acid-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_water      = -5.9_wp     !<  ... water-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28_wp    !<  ... ethanol-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031_wp   !<  ... nonsoluble-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016_wp  !<  ... humus
    !   parameters for temperature dependence
    REAL(wp), PARAMETER   :: temp_p1                    =  0.095_wp   !< first  order temperature dependence [K-1]
    REAL(wp), PARAMETER   :: temp_p2                    = -0.0014_wp  !< second order temperature dependence [K-1]
    !   parameter for precipitation dependence
    REAL(wp), PARAMETER   :: precip_p1                  = -1.21_wp    !< first  order precipitation dependence [yr-1]
    !   parameters for size dependence
    REAL(wp), PARAMETER   :: size_p1                    = -1.71_wp    !< first  order size dependence [cm-1]
    REAL(wp), PARAMETER   :: size_p2                    =  0.86_wp    !< second order size dependence [cm-1]
    REAL(wp), PARAMETER   :: size_p3                    = -0.306_wp   !< size dependence power [ ]

    ! Fractions to spit the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools,
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(wp)   :: A_2_W
    REAL(wp)   :: A_2_E
    REAL(wp)   :: A_2_N
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(wp)   :: W_2_A
    REAL(wp)   :: W_2_E
    REAL(wp)   :: W_2_N
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(wp)   :: E_2_A
    REAL(wp)   :: E_2_W
    REAL(wp)   :: E_2_N
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(wp)   :: N_2_A
    REAL(wp)   :: N_2_W
    REAL(wp)   :: N_2_E
    ! Alternative pathway for N pool which enters the humus pool [ ]
    REAL(wp)   :: N_2_H
    REAL(wp)   :: AWE_2_H

    ! Yasso pools [mol(C) m-2]
    !   Aboveground
    REAL(wp)     :: c_acid_ag
    REAL(wp)     :: c_water_ag
    REAL(wp)     :: c_ethanol_ag
    REAL(wp)     :: c_nonsoluble_ag
    !   Belowground
    REAL(wp)     :: c_acid_bg
    REAL(wp)     :: c_water_bg
    REAL(wp)     :: c_ethanol_bg
    REAL(wp)     :: c_nonsoluble_bg
    REAL(wp)     :: c_humus

    ! Carbon fluxes [mol(C) m-2 yr-1] - Unit indicated by names CAflx (Carbon annual flux)
    !   Aboveground
    REAL(wp)     :: Cflx_yr_from_acid_ag              !< Loss flux of acid soluble pool
    REAL(wp)     :: Cflx_yr_from_water_ag             !< Loss flux of water soluble pool
    REAL(wp)     :: Cflx_yr_from_ethanol_ag           !< Loss flux of ethanol soluble pool
    REAL(wp)     :: Cflx_yr_from_nonsoluble_ag        !< Loss flux of non soluble pool

    REAL(wp)     :: Cflx_yr_2_acid_ag                 !< Gain flux of acid soluble pool
    REAL(wp)     :: Cflx_yr_2_water_ag                !< Gain flux of water soluble pool
    REAL(wp)     :: Cflx_yr_2_ethanol_ag              !< Gain flux of ethanol soluble pool
    REAL(wp)     :: Cflx_yr_2_nonsoluble_ag           !< Gain flux of non soluble pool

    !   Belowground
    REAL(wp)     :: Cflx_yr_from_acid_bg              !< Loss flux of acid soluble pool
    REAL(wp)     :: Cflx_yr_from_water_bg             !< Loss flux of water soluble pool
    REAL(wp)     :: Cflx_yr_from_ethanol_bg           !< Loss flux of ethanol soluble pool
    REAL(wp)     :: Cflx_yr_from_nonsoluble_bg        !< Loss flux of non soluble pool
    REAL(wp)     :: Cflx_yr_from_humus                !< Loss flux of humus pool

    REAL(wp)     :: Cflx_yr_2_acid_bg                 !< Gain flux of acid soluble pool
    REAL(wp)     :: Cflx_yr_2_water_bg                !< Gain flux of water soluble pool
    REAL(wp)     :: Cflx_yr_2_ethanol_bg              !< Gain flux of ethanol soluble pool
    REAL(wp)     :: Cflx_yr_2_nonsoluble_bg           !< Gain flux of non soluble pool

    !   Humus
    REAL(wp)     :: Cflx_yr_2_humus                   !< Gain flux of humus pool
    REAL(wp)     :: Cflx_yr_2_humusAG                 !< Gain flux of humus pool; from ag only
    REAL(wp)     :: Cflx_yr_2_humusBG                 !< Gain flux of humus pool; from bg only
    REAL(wp)     :: Cflx_yr_litter_2_humus            !< Litter influx to humus soluble pool

    ! Respiration is output of Yasso
    REAL(wp)     :: Cflx_yr_respiration               !< Total flux to the atmosphere
    REAL(wp)     :: Cflx_yr_resp_from_litter_ag       !< Flux to the atmosphere; from AWEN ag only
    REAL(wp)     :: Cflx_yr_resp_from_litter_bg       !< Flux to the atmosphere; from AWEN bg only

    ! Decomposition rates [yr-1]
    REAL(wp)     :: d_acid
    REAL(wp)     :: d_water
    REAL(wp)     :: d_ethanol
    REAL(wp)     :: d_nonsoluble
    REAL(wp)     :: d_humus
    !   For methane relevant processes
    REAL(wp)     :: d_anoxic_ag             !< decomposition rate modifier for ag pools, if conditions are anoxic
    REAL(wp)     :: d_anoxic_bg             !< decomposition rate modifier for bg pools, if conditions are anoxic
    REAL(wp)     :: d_peat_acrotelm         !< decomposition rate modifier for peatland C pools acrotelm
    REAL(wp)     :: d_peat_catotelm         !< decomposition rate modifier for peatland C pools catotelm

    REAL(wp)     :: d_temp                  !< term which accounts for the temperature influence on litter decomposition
    REAL(wp)     :: d_precip                !< term which accounts for the precipitation influence on litter decomposition
    REAL(wp)     :: d_size                  !< term which accounts for the litter size influence on litter decomposition

    ! Decomposition rate for N litter green (needed for N cycle)
    REAL(wp)     :: d_litter_green          !< the rate with which the sum of AWEN pools decays
    REAL(wp)     :: Pseudo_litter_green     !< Sum of all AWEN pools

    ! climatic drivers
    REAL(wp)     :: precip                  !< runnning mean of precipitation      [m yr-1]
    REAL(wp)     :: temp                    !< runnning mean of 2m air temperature [C]

    REAL(wp)     :: redFactor               !< dummy for redFact_Nlimit

    ! anoxic & peat decomposition modifier
    LOGICAL      :: anoxic_flag, peat_flag

    ! time stepping: yasso is parameterized based on annual time
    ! steps, thus the fluxes have to be scaled up to annual rates.
    REAL(wp)     :: years_per_ts, ts_per_year

    ts_per_year = REAL(sec_per_year,wp) / delta_time_C
    years_per_ts = 1.0_wp/ts_per_year

    ! Treatment of optional parameters
    !   Set exchange coefficients
    IF (PRESENT(AWEN2H)) THEN
      AWE_2_H = AWEN2H
      N_2_H   = AWEN2H
    END IF
    ! Default values in [yr-1]
    AWE_2_H  = 0.0045_wp; IF (PRESENT(AWE2H))  AWE_2_H  = AWE2H
    N_2_H    = 0.0045_wp; IF (PRESENT(  N2H))    N_2_H  =   N2H
    A_2_W    = 0.99_wp;   IF (PRESENT(  A2W))    A_2_W  =   A2W
    A_2_E    = 0.00_wp;   IF (PRESENT(  A2E))    A_2_E  =   A2E
    A_2_N    = 0.00_wp;   IF (PRESENT(  A2N))    A_2_N  =   A2N
    W_2_A    = 0.48_wp;   IF (PRESENT(  W2A))    W_2_A  =   W2A
    W_2_E    = 0.00_wp;   IF (PRESENT(  W2E))    W_2_E  =   W2E
    W_2_N    = 0.015_wp;  IF (PRESENT(  W2N))    W_2_N  =   W2N
    E_2_A    = 0.01_wp;   IF (PRESENT(  E2A))    E_2_A  =   E2A
    E_2_W    = 0.00_wp;   IF (PRESENT(  E2W))    E_2_W  =   E2W
    E_2_N    = 0.95_wp;   IF (PRESENT(  E2N))    E_2_N  =   E2N
    N_2_A    = 0.83_wp;   IF (PRESENT(  N2A))    N_2_A  =   N2A
    N_2_W    = 0.01_wp;   IF (PRESENT(  N2W))    N_2_W  =   N2W
    N_2_E    = 0.02_wp;   IF (PRESENT(  N2E))    N_2_E  =   N2E

    anoxic_flag = .FALSE.; IF (PRESENT(is_anoxic)) anoxic_flag = is_anoxic
    peat_flag   = .FALSE.; IF (PRESENT(is_peat))   peat_flag   = is_peat
    redFactor   =  1.0_wp
    IF (PRESENT(redFact_Nlimit)) THEN
      ! nitrogen cycle active & we know already the Nlimitation factor
      redFactor = redFact_Nlimit
    END IF

    ! set up anoxic decomposition rate modifiers
    d_anoxic_ag = 1._wp
    d_anoxic_bg = 1._wp
    IF (anoxic_flag) THEN
      ! anoxic conditions are present, reduce decomposition rates of AWEN pools
      ! following Wania et al. 2009, https://doi.org/10.1029/2008GB003413)
      d_anoxic_ag = 0.35_wp
      d_anoxic_bg = 0.35_wp
    END IF

    ! set up peat decomposition rate modifiers
    d_peat_acrotelm = 1._wp
    d_peat_catotelm = 1._wp
    IF (peat_flag) THEN
      d_peat_acrotelm = 1._wp                   ! no decomposition in litter pools modifier -- we reduce GPP instead
      ! decomposition timescale in LPJ is 7.e-5 -> modifier is applied to (ref_decomp_rate_humus * d_anoxic_bg)
      d_peat_catotelm = 0.125_wp                ! decomposition in catotelm is reduced  -> 0.0016*0.35*0.125 = 7e-5
    END IF

    ! 0. Preparations

    !   Change units of climatic forcing variables
    precip = Weather(2)*sec_per_year/1000._wp  ! mm/s -> m/a
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
    Pseudo_litter_green = c_acid_ag          &
                        + c_water_ag         &
                        + c_ethanol_ag       &
                        + c_nonsoluble_ag    &
                        + c_acid_bg          &
                        + c_water_bg         &
                        + c_ethanol_bg       &
                        + c_nonsoluble_bg


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 1. Calculate decomposition rates

    !   Litter size dependence of decomposition -- no effect if WoodlitterSize = 0.0
    d_size   = MIN(1.0_wp,(1.0_wp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_wp)**size_p3)
    !   Temperature dependence of decomposition
    d_temp   = EXP(temp_p1*temp + temp_p2*temp**2.0_wp)
    !   Precipitation dependence of decomposition (not for peat)
    IF (peat_flag) THEN
      d_precip = 1.0_wp
    ELSE
      d_precip = 1.0_wp - EXP(precip_p1*precip)
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
                                  + Cflx_yr_from_ethanol_ag    * E_2_A  &
                                  + Cflx_yr_from_nonsoluble_ag * N_2_A)

    Cflx_yr_2_water_ag        = ABS(Cflx_yr_from_acid_ag       * A_2_W  &
                                  + Cflx_yr_from_ethanol_ag    * E_2_W  &
                                  + Cflx_yr_from_nonsoluble_ag * N_2_W)

    Cflx_yr_2_ethanol_ag      = ABS(Cflx_yr_from_acid_ag       * A_2_E  &
                                  + Cflx_yr_from_water_ag      * W_2_E  &
                                  + Cflx_yr_from_nonsoluble_ag * N_2_E)

    Cflx_yr_2_nonsoluble_ag   = ABS(Cflx_yr_from_acid_ag       * A_2_N  &
                                  + Cflx_yr_from_water_ag      * W_2_N  &
                                  + Cflx_yr_from_ethanol_ag    * E_2_N)

    Cflx_yr_2_acid_bg         = ABS(Cflx_yr_from_water_bg      * W_2_A  &
                                  + Cflx_yr_from_ethanol_bg    * E_2_A  &
                                  + Cflx_yr_from_nonsoluble_bg * N_2_A)

    Cflx_yr_2_water_bg        = ABS(Cflx_yr_from_acid_bg       * A_2_W  &
                                  + Cflx_yr_from_ethanol_bg    * E_2_W  &
                                  + Cflx_yr_from_nonsoluble_bg * N_2_W)

    Cflx_yr_2_ethanol_bg      = ABS(Cflx_yr_from_acid_bg       * A_2_E  &
                                  + Cflx_yr_from_water_bg      * W_2_E  &
                                  + Cflx_yr_from_nonsoluble_bg * N_2_E)

    Cflx_yr_2_nonsoluble_bg   = ABS(Cflx_yr_from_acid_bg       * A_2_N  &
                                  + Cflx_yr_from_water_bg      * W_2_N  &
                                  + Cflx_yr_from_ethanol_bg    * E_2_N)

    Cflx_yr_2_humusAG         = ABS(Cflx_yr_from_acid_ag                &
                                  + Cflx_yr_from_water_ag               &
                                  + Cflx_yr_from_ethanol_ag             &
                                  ) * AWE_2_H                           &
                                  +                                     &
                                ABS(Cflx_yr_from_nonsoluble_ag)         &
                                    * N_2_H

    Cflx_yr_2_humusBG         = ABS(Cflx_yr_from_acid_bg                &
                                  + Cflx_yr_from_water_bg               &
                                  + Cflx_yr_from_ethanol_bg             &
                                    ) * AWE_2_H                         &
                                    +                                   &
                                ABS(Cflx_yr_from_nonsoluble_bg)         &
                                    * N_2_H

    Cflx_yr_2_humus           = Cflx_yr_2_humusAG + Cflx_yr_2_humusBG

    !   the remaining fractions of the loss fluxes enter the atmosphere as respiration
!    Cflx_yr_respiration        =  (Cflx_yr_from_acid_ag + Cflx_yr_from_acid_bg)             &
!                              * (1.0_wp - A_2_W - A_2_E - A_2_N - AWE_2_H)            &
!                              + (Cflx_yr_from_water_ag + Cflx_yr_from_water_bg)           &
!                              * (1.0_wp - W_2_A - W_2_E - W_2_N - AWE_2_H)            &
!                              + (Cflx_yr_from_ethanol_ag + Cflx_yr_from_ethanol_bg)       &
!                              * (1.0_wp - E_2_A - E_2_W - E_2_N - AWE_2_H)            &
!                              + (Cflx_yr_from_nonsoluble_ag + Cflx_yr_from_nonsoluble_bg) &
!                              * (1.0_wp - N_2_A - N_2_W - N_2_E - N_2_H)              &
!                              + Cflx_yr_from_humus

    !   litter ag & bg respiration (needed for N cycle)
    Cflx_yr_resp_from_litter_ag = Cflx_yr_from_acid_ag       * (1.0_wp - A_2_W - A_2_E - A_2_N - AWE_2_H) &
                                + Cflx_yr_from_water_ag      * (1.0_wp - W_2_A - W_2_E - W_2_N - AWE_2_H) &
                                + Cflx_yr_from_ethanol_ag    * (1.0_wp - E_2_A - E_2_W - E_2_N - AWE_2_H) &
                                + Cflx_yr_from_nonsoluble_ag * (1.0_wp - N_2_A - N_2_W - N_2_E -   N_2_H)

    Cflx_yr_resp_from_litter_bg = Cflx_yr_from_acid_bg       * (1.0_wp - A_2_W - A_2_E - A_2_N - AWE_2_H) &
                                + Cflx_yr_from_water_bg      * (1.0_wp - W_2_A - W_2_E - W_2_N - AWE_2_H) &
                                + Cflx_yr_from_ethanol_bg    * (1.0_wp - E_2_A - E_2_W - E_2_N - AWE_2_H) &
                                + Cflx_yr_from_nonsoluble_bg * (1.0_wp - N_2_A - N_2_W - N_2_E -   N_2_H)

    Cflx_yr_respiration         = Cflx_yr_resp_from_litter_ag                      &
                                + Cflx_yr_resp_from_litter_bg                      &
                                + Cflx_yr_from_humus

    ! 3. update Yasso pools

    !   3a. redistribution among yasso pools
    c_acid_ag        = MAX(0.0_wp,                                         &
                          c_acid_ag                                        & ! old pool
                          + (Cflx_yr_from_acid_ag                          & ! incoming flux from AWEN pools
                            + Cflx_yr_2_acid_ag) * years_per_ts)
    c_water_ag       = MAX(0.0_wp,                                         & ! and so on .....
                          c_water_ag                                       &
                          + (Cflx_yr_from_water_ag                         &
                            + Cflx_yr_2_water_ag) * years_per_ts)
    c_ethanol_ag     = MAX(0.0_wp,                                         &
                          c_ethanol_ag                                     &
                          + (Cflx_yr_from_ethanol_ag                       &
                            + Cflx_yr_2_ethanol_ag) * years_per_ts)
    c_nonsoluble_ag  = MAX(0.0_wp,                                         &
                          c_nonsoluble_ag                                  &
                          + (Cflx_yr_from_nonsoluble_ag                    &
                            + Cflx_yr_2_nonsoluble_ag) * years_per_ts)

    c_acid_bg        = MAX(0.0_wp,                                         &
                          c_acid_bg                                        & ! old pool
                          + (Cflx_yr_from_acid_bg                          & ! incoming flux from AWEN pools
                            + Cflx_yr_2_acid_bg) * years_per_ts)
    c_water_bg       = MAX(0.0_wp,                                         & ! ...
                          c_water_bg                                       &
                          + (Cflx_yr_from_water_bg                         &
                            + Cflx_yr_2_water_bg) * years_per_ts           &
                          + c_2_rootExudates)                               ! exudates are carbohydrates only and belowground
    c_ethanol_bg     = MAX(0.0_wp,                                         &
                          c_ethanol_bg                                     &
                          + (Cflx_yr_from_ethanol_bg                       &
                            + Cflx_yr_2_ethanol_bg) * years_per_ts)
    c_nonsoluble_bg  = MAX(0.0_wp,                                         &
                          c_nonsoluble_bg                                  &
                          + (Cflx_yr_from_nonsoluble_bg                    &
                            + Cflx_yr_2_nonsoluble_bg) * years_per_ts)

    c_humus          = MAX(0.0_wp,                                         &
                             c_humus                                       &
                             + (Cflx_yr_from_humus                         &
                               + Cflx_yr_2_humus) * years_per_ts)

    !   3b. distribute litter
    !    ag
    CALL distribute_yasso_litter(c_acid_ag, c_water_ag, c_ethanol_ag, c_nonsoluble_ag, c_humus, &
  &                          c_2_litter,          fract_aboveground, Lit_coefV)
    !    bg
    CALL distribute_yasso_litter(c_acid_bg, c_water_bg, c_ethanol_bg, c_nonsoluble_bg, c_humus, &
  &                          c_2_litter, 1.0_wp - fract_aboveground, Lit_coefV)

    ! Cflx_litter_2_humus is not only required for distribution, but possibly also as output for the eq script
    Cflx_yr_litter_2_humus = Lit_coefV(LCT_HUMUS) * c_2_litter * ts_per_year

    ! compute d_litter_green; this is given by the change in AWEN pools
    !          Flx                 Pool(t)                   Pool(t+1)
    !    d_P = ----     &   Flx = ----------   --->   d_P = -----------  - 1
    !          Pool(t)             Pool(t+1)                 Pool(t)

    IF (Pseudo_litter_green > 0.0_wp) THEN
       d_litter_green = (  c_acid_ag          &
                         + c_water_ag         &
                         + c_ethanol_ag       &
                         + c_nonsoluble_ag    &
                         + c_acid_bg          &
                         + c_water_bg         &
                         + c_ethanol_bg       &
                         + c_nonsoluble_bg)   &
                         / Pseudo_litter_green - 1._wp
    ELSE
      d_litter_green = 0.0_wp
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
  &                                litter, fraction, coefficients)

    !$ACC ROUTINE SEQ

    USE mo_carbon_constants, ONLY:   LCT_ACID, LCT_WATER, LCT_ETHANOL, LCT_NONSOLUBLE, LCT_HUMUS

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(INOUT)             :: c_acid        !< yasso acid pool
    REAL(wp), INTENT(INOUT)             :: c_water       !< yasso water pool
    REAL(wp), INTENT(INOUT)             :: c_ethanol     !< yasso ethanol pool
    REAL(wp), INTENT(INOUT)             :: c_nonsoluble  !< yasso nonsoluble pool
    REAL(wp), INTENT(INOUT)             :: c_humus       !< yasso humus pool
    REAL(wp), INTENT(IN)                :: litter        !< litter input (above and belowground) [mol(C) m-2 day-1]
    REAL(wp), INTENT(IN)                :: fraction      !< fraction of litter put into yasso pools []
    REAL(wp), DIMENSION(5), INTENT(IN)  :: coefficients  !< litter coefficient (acid vs water vs ...) [ ]
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
    REAL(wp), INTENT(INOUT) :: this_pool     !< target yasso pool
    REAL(wp), INTENT(IN)    :: litter        !< litter input  (above and belowground) [mol(C) m-2 day-1]
    REAL(wp), INTENT(IN)    :: fraction      !< fraction of litter put into this pool (BG vs AG) [ ]
    REAL(wp), INTENT(IN)    :: coefficient   !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !
    this_pool           = MAX(0.0_wp, this_pool + (fraction * coefficient * litter))

  END SUBROUTINE add_litter_to_yasso_pool

  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2
  !
  SUBROUTINE get_per_tile(  &
  & Cpool_per_tile,                 & ! out
  & Cpool_per_canopy,               & ! in
  & veg_fract_correction,           & ! in
  & fract_fpc_max)                    ! in

    REAL(wp), INTENT(IN), DIMENSION(:) :: &
  & Cpool_per_canopy, &
  & fract_fpc_max,    &
  & veg_fract_correction
    REAL(wp), INTENT(OUT), DIMENSION(:) :: &
  & Cpool_per_tile

    INTEGER :: ic

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, SIZE(fract_fpc_max)
      Cpool_per_tile(ic) = Cpool_per_canopy(ic) * veg_fract_correction(ic) * fract_fpc_max(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE get_per_tile

#endif
END MODULE mo_carbon_shared
