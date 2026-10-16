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
  USE mo_exception,        ONLY: finish
  USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround,    &
    & Y_ACID, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS, ABOVE_GROUND, BELOW_GROUND, &
    & molarMassCO2_kg

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_plant_c, calc_soil_c, relocate_carbon_fire, relocate_carbon_damage, &
    & get_per_tile, add_litter_to_yasso_pool, sum_up_carbon_layers

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_process'

CONTAINS

  ! Calculations to allocate NPP to the plant C pools as well as various respiration and litter pools.
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_plant_c( &
    ! Inputs
    & delta_time_C,               &
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
    & pheno_type,                 &
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
    !

    !! Names of carbon variables and arguments are:
    !!   1) Cflx_*  : Carbon flux [mol(C) m-2 s-1]
    !!   2) NPP*    : as 1)
    !!   3) c_*_2_* : Carbon transport [mol(C) m-2 carbon-time-step-1], where carbon-time-step is normally 1 day
    !!   4) c_*     : (if not covered by 3)) Carbon pool [mol(C) m-2]
    !! All [m-2] are refere to canopy area

    ! Input arguments
    REAL(lp), INTENT(in)    :: delta_time_C            !< Carbon time step [s]
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
    INTEGER,  INTENT(in)    :: pheno_type              !< Phenology type (checks for crop == 5)
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

    LOGICAL  :: is_crop                 ! logical crop mask

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

    ! Preparations
    !-------------------------------------------

    IF (pheno_type == 5) THEN
       is_crop = .TRUE.
    ELSE
       is_crop = .FALSE.
    END IF

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

  !
  !> This routine handles the impact of fire on the carbon pools.
  !>
  !> We assume that the above ground litter carbon (yasso ag1 and ag2 pools), is completely released
  !> to the atmosphere while the carbon of the living plants (c_green, c_reserve, c_woods) is partly
  !> released to the atmosphere and partly relocated into the litter pools.
  !
  SUBROUTINE relocate_carbon_fire(       &
    & delta_time_C,                      & ! in
    & burned_fract,                      & ! in
    & fire_fract_wood_2_atmos,           & ! in
    & LeafLit_coef,                      & ! in
    & WoodLit_coef,                      & ! in
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
    REAL(lp), INTENT(in)    :: burned_fract(:)              !< Fraction of the vegetated area of each tile burned till the last
                                                            !< call of this routine
    REAL(lp), INTENT(in)    :: fire_fract_wood_2_atmos      !< Fraction of above ground wood immediately emitted to atmosphere
                                                            !< in case of fire
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)              !< Factor to spread non woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)              !< Factor to spread woody litterfall into yasso pools [ ]
    REAL(lp), INTENT(inout) :: c_green(:)                   !< Amount of green carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_woods(:)                   !< Amount of wood carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_reserve(:)                 !< Amount of reserve carbon [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: cflux_dist_green_2_soil(:)   !< Carbon relocated by disturbances from non-woody biomass
                                                            !< to the green litter pools [mol(C)/m^2(canopy)s-1]
    REAL(lp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   !< Carbon relocated by disturbances from woody biomass to
                                                            !< the wood litter pools [mol(C)/m^2(canopy)s-1]
    !   YASSO
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    REAL(lp), INTENT(inout) :: c_acid_ag1(:)                !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_ag1(:)               !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_ag1(:)             !< Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag1(:)          !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp), INTENT(inout) :: c_acid_bg1(:)                !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg1(:)               !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg1(:)             !< Ethanol soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg1(:)          !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_1(:)                 !< Humus pool from green litter [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    REAL(lp), INTENT(inout) :: c_acid_ag2(:)                !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_ag2(:)               !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_ag2(:)             !< Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_ag2(:)          !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp), INTENT(inout) :: c_acid_bg2(:)                !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_water_bg2(:)               !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_ethanol_bg2(:)             !< Ethanol soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_nonsoluble_bg2(:)          !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(inout) :: c_humus_2(:)                 !< Humus pool from woody litter [mol(C)/m^2(canopy)]
    REAL(lp), INTENT(out)   :: co2flux_fire_all_2_atm(:)    !< CO2 released by fire [kg(CO2)m-2(canopy)s-1]

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
      cflux_dist_woods_2_soil(ic) = cflux_dist_woods_2_soil(ic) + &
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
      c_acid_ag1(ic)       = fraction * c_acid_ag1(ic)
      c_water_ag1(ic)      = fraction * c_water_ag1(ic)
      c_ethanol_ag1(ic)    = fraction * c_ethanol_ag1(ic)
      c_nonsoluble_ag1(ic) = fraction * c_nonsoluble_ag1(ic)

      c_acid_ag2(ic)       = fraction * c_acid_ag2(ic)
      c_water_ag2(ic)      = fraction * c_water_ag2(ic)
      c_ethanol_ag2(ic)    = fraction * c_ethanol_ag2(ic)
      c_nonsoluble_ag2(ic) = fraction * c_nonsoluble_ag2(ic)

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

  !
  !> This routine handles the impact of vegetation damage on the carbon pools.
  !>
  !< In this routine we model the impact of vegetation damage (i.e. wind through) on the carbon pools:
  !< Depending on the disturbed area fraction carbon from the living vegetation pools (c_green,
  !< c_reserve and c_woods) is relocated  to the litter and soil pools.
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
    REAL(lp),INTENT(inout) :: c_acid_ag1(:)                !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_ag1(:)               !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_ag1(:)             !< Ethanol soluble green litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_ag1(:)          !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp),INTENT(inout) :: c_acid_bg1(:)                !< Acid soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_bg1(:)               !< Water soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_bg1(:)             !< Ethanol soluble green litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_bg1(:)          !< Non-soluble green litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_humus_1(:)                 !< Humus pool from green litter [mol(C)/m^2(canopy)]
    ! Size class 2: wood
    !   above ground C pools
    REAL(lp),INTENT(inout) :: c_acid_ag2(:)                !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_ag2(:)               !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_ag2(:)             !< Ethanol soluble woody litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_ag2(:)          !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(lp),INTENT(inout) :: c_acid_bg2(:)                !< Acid soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_water_bg2(:)               !< Water soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_ethanol_bg2(:)             !< Ethanol soluble woody litter pool [mol(C)/m2(canopy)]
    REAL(lp),INTENT(inout) :: c_nonsoluble_bg2(:)          !< Non-soluble woody litter pool [mol(C)/m^2(canopy)]
    REAL(lp),INTENT(inout) :: c_humus_2(:)                 !< Humus pool from woody litter [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(lp)  :: non_woody_litter, fraction  ! Helper
    INTEGER   :: ic, nc

    ! Initializations
    nc = SIZE(damaged_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) &
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
      CALL distribute_yasso_litter(c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), &
        &                          c_nonsoluble_ag1(ic), c_humus_1(ic),                &
        &                          non_woody_litter, fraction, LeafLit_coef)
      !  Below ground
      fraction = (1._lp - fract_green_aboveGround) * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), &
        &                          c_nonsoluble_bg1(ic), c_humus_1(ic),                &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the carbon of damaged woody vegetation to the yasso litter and soil pools
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


  !
  !> Compute soil carbon dynamics.
  !>
  !> Carbon fluxes from litter and root exudates are added to the respective carbon pools,
  !> which are then decomposed into humus according to the Yasso model.
  !
  SUBROUTINE calc_soil_c(vertical_yasso, nsoil, delta_time_C,                               &
    & LeafLit_coef, WoodLit_coef, MiscPar, temp, precip,                                    &
    & Cflx_wood_2_litter, Cflx_green_2_litter, Cflx_2_root_exudates,                        &
    & green_cpools_ag, woody_cpools_ag, Cflx_soil_respiration,                              &
    & c_decomp_humus_1_sum, c_decomp_humus_2_sum, c_into_humus_1_sum, c_into_humus_2_sum,   &
    & t_soil_sl, wtr_soil_rel_sl, n_thaw_cycle, soil_depth, soil_depth_sl, max_depth,       &
    & max_depth_sl, layer_depth, layer_thick, thaw_depth, permafrost_period_sl,             &
    & green_cpools_bg,    green_cpools_humus,    woody_cpools_bg,    woody_cpools_humus,    &
    & green_cpools_bg_sl, green_cpools_humus_sl, woody_cpools_bg_sl, woody_cpools_humus_sl  &
    & )

    !$ACC ROUTINE SEQ

    USE mo_carbon_constants, ONLY: WOODLITSIZE, ROOT_DIST_BETA

    LOGICAL,  INTENT(in)    :: vertical_yasso        !< Yasso scheme switch
    INTEGER,  INTENT(in)    :: nsoil                 !< Number of soil layers in vertical yasso scheme []
    REAL(lp), INTENT(in)    :: delta_time_C          !< Carbon time step [s]
    !   Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(lp), INTENT(in)    :: LeafLit_coef(:)       !< Leaf litter fraction going to specific pools
    !   Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(lp), INTENT(in)    :: WoodLit_coef(:)       !< Woody litter fraction going to specific pools
    REAL(lp), INTENT(in)    :: MiscPar(:)            !< miscellaneous parameters

    ! Carbon inputs [mol m-2 s-1]:
    REAL(lp), INTENT(inout) :: Cflx_wood_2_litter    !< Carbon flux from woody litter to soil [mol(C) m-2 s-1]
    REAL(lp), INTENT(inout) :: Cflx_green_2_litter   !< Carbon flux from green litter to soil [mol(C) m-2 s-1]
    REAL(lp), INTENT(inout) :: Cflx_2_root_exudates  !< Carbon flux from root exudates to soil [mol(C) m-2 s-1]

    ! Input meteorology for Yasso
    REAL(lp), INTENT(in)    :: temp                  !< 15-day mean surface temperature [K]
    REAL(lp), INTENT(in)    :: precip                !< 15-day mean precipitation [kg m-2 s-1]

    ! Soil carbon pools [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout) :: green_cpools_ag(:)    !< Aboveground carbon pools of green litter [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout) :: woody_cpools_ag(:)    !< Aboveground carbon pools of woody litter [mol(C) m-2(canopy)]

    ! Decomposistion fluxes
    REAL(lp), INTENT(inout) :: c_decomp_humus_1_sum  !< Annual sum of carbon of humus (leaf) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_decomp_humus_2_sum  !< Annual sum of carbon of humus (wood) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_into_humus_1_sum    !< Annual sum of carbon decomposed into humus (leaf) [mol(C) m-2 yr-1]
    REAL(lp), INTENT(inout) :: c_into_humus_2_sum    !< Annual sum of carbon decomposed into humus (wood) [mol(C) m-2 yr-1]

    REAL(lp), INTENT(out)   :: Cflx_soil_respiration !< Soil (=heterotrophic) respiration rate  [mol(C) m-2 s-1]
                                                     ! Note: this is a loss, therefore negative!

    ! Optional arguments for either single or multilayer yasso
    REAL(lp), INTENT(in),    OPTIONAL :: t_soil_sl(:)             !< soil layer temperatures [K]
    REAL(lp), INTENT(in),    OPTIONAL :: wtr_soil_rel_sl(:)       !< relative soil moisture []
    REAL(lp), INTENT(in),    OPTIONAL :: n_thaw_cycle             !< number of freeze-thaw cycles in previous year
    REAL(lp), INTENT(in),    OPTIONAL :: soil_depth               !< Total soil depth until bedrock [m]
    REAL(lp), INTENT(in),    OPTIONAL :: soil_depth_sl(:)         !< Soil depth within each layer [m]
    REAL(lp), INTENT(in),    OPTIONAL :: max_depth                !< maximum available depth [m], either soil or thaw
    REAL(lp), INTENT(in),    OPTIONAL :: max_depth_sl(:)          !< maximum available depth [m], either soil or thaw, per layer
    REAL(lp), INTENT(in),    OPTIONAL :: layer_depth(:)           !< Depths of soil layer center [m]
    REAL(lp), INTENT(in),    OPTIONAL :: layer_thick(:)           !< Thickness of soil layer [m]
    REAL(lp), INTENT(in),    OPTIONAL :: thaw_depth               !< thaw depth [m]
    REAL(lp), INTENT(in),    OPTIONAL :: permafrost_period_sl(:)  !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(inout), OPTIONAL :: green_cpools_bg(:)       !< Belowground carbon pools of green litter [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout), OPTIONAL :: green_cpools_humus       !< Carbon pool of humus (leaf) [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout), OPTIONAL :: woody_cpools_bg(:)       !< Belowground carbon pools of woody litter [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout), OPTIONAL :: woody_cpools_humus       !< Carbon pool of humus (wood) [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout), OPTIONAL :: green_cpools_bg_sl(:,:)  !< Belowground carbon pools of green litter [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout), OPTIONAL :: green_cpools_humus_sl(:) !< Carbon pool of humus (leaf) [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout), OPTIONAL :: woody_cpools_bg_sl(:,:)  !< Belowground carbon pools of woody litter [mol(C) m-2(canopy)]
    REAL(lp), INTENT(inout), OPTIONAL :: woody_cpools_humus_sl(:) !< Carbon pool of humus (wood) [mol(C) m-2(canopy)]

    ! local variables for Yasso
    INTEGER  :: is                  !< Soil layer index []
    INTEGER  :: ip                  !< Pool type index []
    INTEGER  :: npools              !< Maximum number of above and below ground AWEN pools []
    REAL(lp) :: Weather(5)          !< Vector containing climate metrics for Yasso
    REAL(lp) :: LeafLitterSize      !< Leaf litter size (always 0 -- > size has no effect on decomposition)
    REAL(lp) :: WoodLitterSize      !< Wood litter size (PFT dependent)
    REAL(lp) :: ZeroFlux            !< Flux that is zero for a specific pool (green,wood,above,below)
    REAL(lp) :: root_dist_sl(nsoil) !< Root distribution factor for each soil layer []
    REAL(lp) :: diff_coeff(nsoil)   !< Vertical carbon flux diffusion coefficient [m2 s-1]
    REAL(lp) :: Cflx_soil_respiration_layer
    REAL(lp) :: Cflx_into_humus_layer
    REAL(lp) :: Cflx_decomp_humus_layer

    ! Preparations
    !-------------------------------------------

    Cflx_soil_respiration = 0._lp
    Weather(:)            = 0._lp
    LeafLitterSize        = 0._lp
    WoodLitterSize        = MiscPar(WOODLITSIZE)
    ZeroFlux              = 0._lp

    Weather(1)            = temp
    Weather(2)            = precip

    npools = SIZE(green_cpools_ag(:))

    IF (vertical_yasso) THEN
      root_dist_sl(:) = determine_root_dist(max_depth_sl, layer_depth, MiscPar(ROOT_DIST_BETA), &
        &                                   t_soil_sl, permafrost_period_sl)
    END IF

    ! ============== COMPUTE SOIL RESPIRATION for single layer yasso ============================================
    !
    IF (.NOT. vertical_yasso) THEN
      ! Compute yasso with single layer below ground pools
      ! Above ground green (leaf) litter
      CALL yasso(vertical_yasso, ABOVE_GROUND, delta_time_C, Weather(:), Cflx_green_2_litter * delta_time_C,    &
        &        LeafLit_coef(:), LeafLitterSize, fract_green_aboveGround,                      &
        &        Cflx_2_root_exudates * delta_time_C, green_cpools_ag(:), green_cpools_humus,   &
        &        Cflx_soil_respiration_layer, Cflx_into_humus_layer, Cflx_decomp_humus_layer)
      c_into_humus_1_sum    = c_into_humus_1_sum    + Cflx_into_humus_layer
      Cflx_soil_respiration = Cflx_soil_respiration + Cflx_soil_respiration_layer / delta_time_C

      ! Above ground woody litter
      CALL yasso(vertical_yasso, ABOVE_GROUND, delta_time_C, Weather(:), Cflx_wood_2_litter * delta_time_C,     &
        &        WoodLit_coef(:), WoodLitterSize, fract_wood_aboveGround,                       &
        &        ZeroFlux, woody_cpools_ag(:), woody_cpools_humus,                              &
        &        Cflx_soil_respiration_layer, Cflx_into_humus_layer, Cflx_decomp_humus_layer)
      c_into_humus_2_sum    = c_into_humus_2_sum    + Cflx_into_humus_layer
      Cflx_soil_respiration = Cflx_soil_respiration + Cflx_soil_respiration_layer / delta_time_C

      ! Below ground green (leaf) litter
      CALL yasso(vertical_yasso, BELOW_GROUND, delta_time_C, Weather(:), Cflx_green_2_litter * delta_time_C,    &
        &        LeafLit_coef(:), LeafLitterSize, fract_green_aboveGround,                      &
        &        Cflx_2_root_exudates * delta_time_C, green_cpools_bg(:), green_cpools_humus,   &
        &        Cflx_soil_respiration_layer, Cflx_into_humus_layer, Cflx_decomp_humus_layer)
      c_into_humus_1_sum    = c_into_humus_1_sum    + Cflx_into_humus_layer
      c_decomp_humus_1_sum  = c_decomp_humus_1_sum  + Cflx_decomp_humus_layer
      Cflx_soil_respiration = Cflx_soil_respiration + Cflx_soil_respiration_layer / delta_time_C

      ! Below ground woody litter
      CALL yasso(vertical_yasso, BELOW_GROUND, delta_time_C, Weather(:), Cflx_wood_2_litter * delta_time_C,     &
        &        WoodLit_coef(:), WoodLitterSize, fract_wood_aboveGround,                       &
        &        ZeroFlux, woody_cpools_bg(:), woody_cpools_humus,                              &
        &        Cflx_soil_respiration_layer, Cflx_into_humus_layer, Cflx_decomp_humus_layer)
      c_into_humus_2_sum    = c_into_humus_2_sum    + Cflx_into_humus_layer
      c_decomp_humus_2_sum  = c_decomp_humus_2_sum  + Cflx_decomp_humus_layer
      Cflx_soil_respiration = Cflx_soil_respiration + Cflx_soil_respiration_layer / delta_time_C

    ! ============== COMPUTE SOIL RESPIRATION for multi-layer yasso =============================================
    !
    ELSE
      ! Compute yasso with multiple soil layers for the below ground pools
      ! Above ground green (leaf) litter
      CALL yasso(vertical_yasso, ABOVE_GROUND, delta_time_C, Weather(:), Cflx_green_2_litter * delta_time_C,    &
        &        LeafLit_coef(:), LeafLitterSize, fract_green_aboveGround,                      &
        &        ZeroFlux, green_cpools_ag(:), green_cpools_humus_sl(1),                        &
        &        Cflx_soil_respiration_layer, Cflx_into_humus_layer, Cflx_decomp_humus_layer)
      c_into_humus_1_sum    = c_into_humus_1_sum    + Cflx_into_humus_layer
      Cflx_soil_respiration = Cflx_soil_respiration + Cflx_soil_respiration_layer / delta_time_C

      ! Above ground woody litter
      CALL yasso(vertical_yasso, ABOVE_GROUND, delta_time_C, Weather(:), Cflx_wood_2_litter * delta_time_C,     &
        &        WoodLit_coef(:), WoodLitterSize, fract_wood_aboveGround,                       &
        &        ZeroFlux, woody_cpools_ag(:), woody_cpools_humus_sl(1),                        &
        &        Cflx_soil_respiration_layer, Cflx_into_humus_layer, Cflx_decomp_humus_layer)
      c_into_humus_2_sum    = c_into_humus_2_sum    + Cflx_into_humus_layer
      Cflx_soil_respiration = Cflx_soil_respiration + Cflx_soil_respiration_layer / delta_time_C

      DO is = 1, nsoil
        ! Set level depending climate forcing
        Weather(3) = t_soil_sl(is)
        Weather(4) = wtr_soil_rel_sl(is)
        Weather(5) = layer_depth(is)

        ! Below ground green (leaf) litter
        CALL yasso(vertical_yasso, BELOW_GROUND, delta_time_C, Weather(:), Cflx_green_2_litter * delta_time_C,               &
          &        LeafLit_coef(:), LeafLitterSize, fract_green_aboveGround,                                 &
          &        Cflx_2_root_exudates * delta_time_C, green_cpools_bg_sl(is,:), green_cpools_humus_sl(is), &
          &        Cflx_soil_respiration_layer, Cflx_into_humus_layer, Cflx_decomp_humus_layer,              &
          &        root_dist_sl(is))
        c_into_humus_1_sum    = c_into_humus_1_sum    + Cflx_into_humus_layer
        c_decomp_humus_1_sum  = c_decomp_humus_1_sum  + Cflx_decomp_humus_layer
        Cflx_soil_respiration = Cflx_soil_respiration + Cflx_soil_respiration_layer / delta_time_C

        ! Below ground woody litter
        CALL yasso(vertical_yasso, BELOW_GROUND, delta_time_C, Weather(:), Cflx_wood_2_litter * delta_time_C,     &
          &        WoodLit_coef(:), WoodLitterSize, fract_wood_aboveGround,                       &
          &        ZeroFlux, woody_cpools_bg_sl(is,:), woody_cpools_humus_sl(is),                 &
          &        Cflx_soil_respiration_layer, Cflx_into_humus_layer, Cflx_decomp_humus_layer,   &
          &        root_dist_sl(is))
        c_into_humus_2_sum    = c_into_humus_2_sum    + Cflx_into_humus_layer
        c_decomp_humus_2_sum  = c_decomp_humus_2_sum  + Cflx_decomp_humus_layer
        Cflx_soil_respiration = Cflx_soil_respiration + Cflx_soil_respiration_layer / delta_time_C
      END DO

      ! ============== VERTICAL MIXING OF CARBON POOLS ==========================================================
      !
      CALL move_ag_2_bg(nsoil, max_depth, layer_thick, permafrost_period_sl,  &
        &               green_cpools_ag(:), woody_cpools_ag(:),               &
        &               green_cpools_bg_sl(:,:), woody_cpools_bg_sl(:,:),     &
        &               green_cpools_humus_sl(:), woody_cpools_humus_sl(:))

      CALL set_diffusivity(soil_depth, soil_depth_sl, thaw_depth, layer_depth, t_soil_sl, &
        &                  wtr_soil_rel_sl, n_thaw_cycle, permafrost_period_sl, diff_coeff)
      ! AWEN Pools
      DO ip = 1, npools
        CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, &
          & green_cpools_bg_sl(:,ip))
        CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, &
          & woody_cpools_bg_sl(:,ip))
      END DO
      ! Humus Pools
      CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, &
        & green_cpools_humus_sl(:))
      CALL vertical_C_transport(max_depth_sl, diff_coeff, delta_time_C, layer_thick, permafrost_period_sl, &
        & woody_cpools_humus_sl(:))

      ! remove any C that might be below bedrock boundary
      DO ip = 1, npools
        CALL remove_C_below_bedrock(soil_depth, layer_thick, green_cpools_bg_sl(:,ip))
        CALL remove_C_below_bedrock(soil_depth, layer_thick, woody_cpools_bg_sl(:,ip))
      END DO
    END IF

  END SUBROUTINE calc_soil_c


  SUBROUTINE yasso (vertical_yasso, yasso_case, delta_time_C, Weather, c_2_litter,       &
    &               Lit_coef, LitterSize, fract_aboveground, c_2_rootExudates,           &
    &               cpools_awen, c_humus,                                                &
    &               Cflx_soil_respiration, Cflx_into_humus, Cflx_decomp_humus,           &
    &               root_dist                                                            &
    &               )
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

    USE mo_jsb_physical_constants, ONLY: tmelt
    USE mo_carbon_constants, ONLY: sec_per_year, ref_decomp_rate_acid, ref_decomp_rate_water, &
      & ref_decomp_rate_ethanol, ref_decomp_rate_nonsoluble, ref_decomp_rate_humus,           &
      & A_2_W, A_2_E, A_2_N, W_2_A, W_2_E, W_2_N, E_2_A, E_2_W, E_2_N, N_2_A, N_2_W, N_2_E,   &
      & AWE_2_H, N_2_H

    IMPLICIT NONE

    LOGICAL,  INTENT(in)  :: vertical_yasso           !< T: use vertical yasso scheme
    INTEGER,  INTENT(in)  :: yasso_case               !< Yasso case: ABOVE_GROUND or BELOW_GROUND
    REAL(lp), INTENT(in)  :: delta_time_C             !< Carbon time step [s]
    REAL(lp), INTENT(in)  :: Weather(:)               !< climatic drivers: air temperature and precipitation, 15-d means
    REAL(lp), INTENT(in)  :: c_2_litter               !< fresh litter inputs  (above and belowground) [mol(C) m-2 d-1]
    REAL(lp), INTENT(in)  :: Lit_coef(5)              !< fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(lp), INTENT(in)  :: LitterSize               !< Size of leaf or woody litter
    REAL(lp), INTENT(in)  :: fract_aboveground        !< parameter to seperate litter influx into above- and belowground parts [ ]
    REAL(lp), INTENT(in)  :: c_2_rootExudates         !< root exudates flux [mol(C) m-2 d-1]

    REAL(lp), INTENT(inout) :: cpools_awen(:)         !< Above ground carbon pools [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus                !< Humus pool [mol(C) m-2]

    REAL(lp), INTENT(out)   :: Cflx_soil_respiration  !< Soil respiration flux [mol(C) m-2 d-1]
    REAL(lp), INTENT(out)   :: Cflx_into_humus        !< Carbon flux into humus pool [mol(C) m-2 d-1]
    REAL(lp), INTENT(out)   :: Cflx_decomp_humus      !< Carbon flux from humus pool [mol(C) m-2 d-1]

    REAL(lp), OPTIONAL, INTENT(in) :: root_dist       ! Root distribution factor for a certain soil layer [ ]

    CHARACTER(len=*), PARAMETER :: routine = modname//':yasso'

    ! Yasso pools [mol(C) m-2]
    REAL(lp)   :: c_acid
    REAL(lp)   :: c_water
    REAL(lp)   :: c_ethanol
    REAL(lp)   :: c_nonsoluble

    ! Carbon loss and gain fluxes [mol(C) m-2 yr-1]
    REAL(lp)   :: Cflx_yr_from_acid                 !< Loss flux of acid soluble pool
    REAL(lp)   :: Cflx_yr_from_water                !< Loss flux of water soluble pool
    REAL(lp)   :: Cflx_yr_from_ethanol              !< Loss flux of ethanol soluble pool
    REAL(lp)   :: Cflx_yr_from_nonsoluble           !< Loss flux of non soluble pool
    REAL(lp)   :: Cflx_yr_from_humus                !< Loss flux of humus pool

    REAL(lp)   :: Cflx_yr_2_acid                    !< Gain flux of acid soluble pool
    REAL(lp)   :: Cflx_yr_2_water                   !< Gain flux of water soluble pool
    REAL(lp)   :: Cflx_yr_2_ethanol                 !< Gain flux of ethanol soluble pool
    REAL(lp)   :: Cflx_yr_2_nonsoluble              !< Gain flux of non soluble pool
    REAL(lp)   :: Cflx_yr_2_humus                   !< Gain flux of humus pool

    REAL(lp)   :: Cflx_yr_litter_2_humus            !< Litter influx to humus soluble pool

    ! Respiration is output of Yasso
    REAL(lp)   :: Cflx_yr_respiration               !< Total flux to the atmosphere
    REAL(lp)   :: Cflx_yr_resp_from_litter          !< Flux to the atmosphere; from AWEN ag only

    REAL(lp)   :: fract_litter            !< litter influx fraction into above- or belowground pools [ ]

    ! Decomposition rates [yr-1]
    REAL(lp)   :: d_acid
    REAL(lp)   :: d_water
    REAL(lp)   :: d_ethanol
    REAL(lp)   :: d_nonsoluble
    REAL(lp)   :: d_humus

    ! Climate dependencies
    REAL(lp)   :: d_temp                  !< term which accounts for the surface or soil temperature influence on litter decomposition
    REAL(lp)   :: d_moist                 !< term which accounts for the precipitation or soil moisture influence on litter decomposition
    REAL(lp)   :: d_size                  !< term which accounts for the litter size influence on litter decomposition
    REAL(lp)   :: d_pool                  !< term which accounds for the pool type (above,below) for humus decomposition
    REAL(lp)   :: d_depth                 !< term which accounts for the soil water influence on litter decomposition

    ! Decomposition rate for N litter green (needed for N cycle)
    REAL(lp)   :: d_litter_green          !< the rate with which the sum of AWEN pools decays
    REAL(lp)   :: Pseudo_litter_green     !< Sum of all AWEN pools

    ! climatic drivers
    REAL(lp)   :: precip                  !< runnning mean of precipitation      [m yr-1]
    REAL(lp)   :: temp                    !< runnning mean of 2m air temperature [C]
    REAL(lp)   :: f_roots                 !< root distribution factor for a certain soil layer

    ! time stepping: yasso is parameterized based on annual time
    ! steps, thus the fluxes have to be scaled up to annual rates.
    REAL(lp)   :: years_per_ts, ts_per_year

    ts_per_year  = REAL(sec_per_year, lp) / delta_time_C
    years_per_ts = 1.0_lp / ts_per_year


    ! Checks:
    IF (yasso_case /= ABOVE_GROUND .AND. yasso_case /= BELOW_GROUND) THEN
      CALL finish(routine, 'ERROR: yasso_case must be either ABOVE_GROUND or BELOW_GROUND')
    END IF

    IF (PRESENT(root_dist)) THEN
      f_roots = root_dist
    ELSE IF (yasso_case == BELOW_GROUND .AND. vertical_yasso) THEN
      CALL finish(routine, 'Error: root distribution factor required by vertical yasso')
    ELSE
      f_roots = 1._lp
    END IF

    ! initialize yasso pools (AWEN + humus)
    c_acid        = cpools_awen(Y_ACID)
    c_water       = cpools_awen(Y_WATER)
    c_ethanol     = cpools_awen(Y_ETHANOL)
    c_nonsoluble  = cpools_awen(Y_NONSOLUBLE)

    ! Effects of above vs below ground pools
    IF (yasso_case == ABOVE_GROUND) THEN
      fract_litter = fract_aboveground
      d_pool       = 0._lp
    ELSE IF (yasso_case == BELOW_GROUND) THEN
      fract_litter = MAX(0._lp, 1.0_lp - fract_aboveground) * f_roots
      d_pool       = 1._lp
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 1. Calculate decomposition rates

    !   Litter size dependence of decomposition -- no effect if LitterSize = 0.0 (e.g. for leafs)
    d_size   = get_modifier_size(LitterSize)
    !   Temperature dependence of decomposition
    d_temp   = get_modifier_temp(Weather(1))
    !   Precipitation dependence of decomposition
    d_moist = get_modifier_precip(Weather(2))
    IF (yasso_case == BELOW_GROUND .AND. vertical_yasso) THEN ! Multi-layer yasso
      ! Compute further dependencies
      d_temp  = get_modifier_temp(Weather(3))
      d_moist = get_modifier_water(Weather(4))
      d_depth = get_modifier_depth(Weather(5))
    ELSE
      d_depth = 1._lp
    END IF


    !   decomposition rates accounting for temperature, precipitation, litter size or pool type, and nutrient limitation
    d_acid       = ref_decomp_rate_acid       * d_temp * d_moist * d_size * d_depth
    d_water      = ref_decomp_rate_water      * d_temp * d_moist * d_size * d_depth
    d_ethanol    = ref_decomp_rate_ethanol    * d_temp * d_moist * d_size * d_depth
    d_nonsoluble = ref_decomp_rate_nonsoluble * d_temp * d_moist * d_size * d_depth
    d_humus      = ref_decomp_rate_humus      * d_temp * d_moist * d_pool * d_depth

    ! 2a. Calculate carbon fluxes

    !   annual loss fluxes (negative values):
    Cflx_yr_from_acid       = c_acid       * d_acid
    Cflx_yr_from_water      = c_water      * d_water
    Cflx_yr_from_ethanol    = c_ethanol    * d_ethanol
    Cflx_yr_from_nonsoluble = c_nonsoluble * d_nonsoluble
    Cflx_yr_from_humus      = c_humus      * d_humus

    !   annual gain fluxes (positive values):
    !   fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why)
    Cflx_yr_2_acid         = ABS(Cflx_yr_from_water        * W_2_A  &
      &                        + Cflx_yr_from_ethanol      * E_2_A  &
      &                        + Cflx_yr_from_nonsoluble   * N_2_A)

    Cflx_yr_2_water        = ABS(Cflx_yr_from_acid         * A_2_W  &
      &                        + Cflx_yr_from_ethanol      * E_2_W  &
      &                        + Cflx_yr_from_nonsoluble   * N_2_W)

    Cflx_yr_2_ethanol      = ABS(Cflx_yr_from_acid         * A_2_E  &
      &                        + Cflx_yr_from_water        * W_2_E  &
      &                        + Cflx_yr_from_nonsoluble   * N_2_E)

    Cflx_yr_2_nonsoluble   = ABS(Cflx_yr_from_acid         * A_2_N  &
      &                        + Cflx_yr_from_water        * W_2_N  &
      &                        + Cflx_yr_from_ethanol      * E_2_N)

    Cflx_yr_2_humus        = ABS(Cflx_yr_from_acid                  &
      &                        + Cflx_yr_from_water                 &
      &                        + Cflx_yr_from_ethanol               &
      &                          ) * AWE_2_H                        &
      &                    +                                        &
      &                      ABS(Cflx_yr_from_nonsoluble)           &
      &                            * N_2_H

    !   litter respiration
    Cflx_yr_resp_from_litter =                                                 &
      &   Cflx_yr_from_acid       * (1.0_lp - A_2_W - A_2_E - A_2_N - AWE_2_H) &
      & + Cflx_yr_from_water      * (1.0_lp - W_2_A - W_2_E - W_2_N - AWE_2_H) &
      & + Cflx_yr_from_ethanol    * (1.0_lp - E_2_A - E_2_W - E_2_N - AWE_2_H) &
      & + Cflx_yr_from_nonsoluble * (1.0_lp - N_2_A - N_2_W - N_2_E -   N_2_H)

    ! 2b. update Yasso pools

    !   redistribution among yasso pools
    c_acid        = MAX(0.0_lp,                                      &
      &                 c_acid                                       & ! old pool
      &                 + (Cflx_yr_from_acid                         & ! incoming flux from AWEN pools
      &                 + Cflx_yr_2_acid) * years_per_ts)
    c_ethanol     = MAX(0.0_lp,                                      & ! and so on .....
      &                 c_ethanol                                    &
      &                 + (Cflx_yr_from_ethanol                      &
      &                 + Cflx_yr_2_ethanol) * years_per_ts)
    c_nonsoluble  = MAX(0.0_lp,                                      &
      &                 c_nonsoluble                                 &
      &                 + (Cflx_yr_from_nonsoluble                   &
      &                 + Cflx_yr_2_nonsoluble) * years_per_ts)
    IF (yasso_case == ABOVE_GROUND) THEN
      ! Above ground (no humus pool decomposition or root exudates)
      c_water     = MAX(0.0_lp,                                      &
        &               c_water                                      &
        &               + (Cflx_yr_from_water                        &
        &               + Cflx_yr_2_water) * years_per_ts)
      c_humus     = MAX(0.0_lp,                                      &
        &               c_humus                                      &
        &               + Cflx_yr_2_humus * years_per_ts)
    ELSE IF (yasso_case == BELOW_GROUND) THEN
      ! Below ground (with humus pool decomposition and root exudates)
      c_water     = MAX(0.0_lp,                                      &
        &               c_water                                      &
        &               + (Cflx_yr_from_water                        &
        &               + Cflx_yr_2_water) * years_per_ts            &
        &               + c_2_rootExudates * f_roots)                  ! exudates are carbohydrates only
      c_humus     = MAX(0.0_lp,                                      &
        &               c_humus                                      &
        &               + (Cflx_yr_from_humus                        &
        &               + Cflx_yr_2_humus) * years_per_ts)
    END IF

    ! 4.0 Distribute litter to the pools
    CALL distribute_yasso_litter(c_acid, c_water, c_ethanol, c_nonsoluble, c_humus, &
      &                          c_2_litter, fract_litter, Lit_coef)

    Cflx_yr_litter_2_humus = Lit_coef(Y_HUMUS) * c_2_litter * ts_per_year

    ! 4.1 Collect the updated pools into the output arrays
    cpools_awen(Y_ACID)       = c_acid
    cpools_awen(Y_WATER)      = c_water
    cpools_awen(Y_ETHANOL)    = c_ethanol
    cpools_awen(Y_NONSOLUBLE) = c_nonsoluble

    ! 4.2 Compute the total soil respiration flux (from AWEN and humus pools)
    Cflx_yr_respiration   = Cflx_yr_resp_from_litter + Cflx_yr_from_humus
    Cflx_soil_respiration = Cflx_yr_respiration * years_per_ts ! convert back to daily

    ! 4.3 Write diagnostics to determine equilibrium of humus pools
    Cflx_into_humus   = (Cflx_yr_2_humus + Cflx_yr_litter_2_humus) * years_per_ts
    Cflx_decomp_humus = -Cflx_yr_from_humus * years_per_ts

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

    CALL add_litter_to_yasso_pool(c_acid,       litter, fraction, coefficients(Y_ACID))
    CALL add_litter_to_yasso_pool(c_water,      litter, fraction, coefficients(Y_WATER))
    CALL add_litter_to_yasso_pool(c_ethanol,    litter, fraction, coefficients(Y_ETHANOL))
    CALL add_litter_to_yasso_pool(c_nonsoluble, litter, fraction, coefficients(Y_NONSOLUBLE))
    CALL add_litter_to_yasso_pool(c_humus,      litter, fraction, coefficients(Y_HUMUS))

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
#ifndef _OPENACC
  ELEMENTAL &
#endif
  FUNCTION get_modifier_temp(temp_soil, use_Q10) RESULT(d_temp)
    USE mo_carbon_constants,       ONLY: Q10, temp_p1, temp_p2
    USE mo_jsb_physical_constants, ONLY: tmelt

    REAL(lp), INTENT(in)          :: temp_soil
    LOGICAL, OPTIONAL, INTENT(in) :: use_Q10

    LOGICAL :: do_Q10
    REAL(lp):: d_temp

    IF (PRESENT(use_Q10)) THEN
      do_Q10 = use_Q10
    ELSE
      do_Q10 = .FALSE.
    END IF

    IF (.NOT. do_Q10) THEN       ! use YASSO formulation
      d_temp  = EXP(temp_p1*(temp_soil - tmelt) + &                 ! YASSO formulation for temperature dependence
        &           temp_p2*(temp_soil - tmelt)**2.0_lp)
    ELSE                        ! use alternative Q10 formulation
      d_temp  = Q10**((temp_soil - (tmelt + 10._lp))/10._lp)        ! Temperature modifier as CLM
    END IF
  END FUNCTION get_modifier_temp

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_water(water_soil) RESULT(d_water)
    USE mo_carbon_constants,  ONLY: moist_p_a, moist_p_b, moist_p_d

    REAL(lp), INTENT(in)  :: water_soil
    REAL(lp)              :: d_water

    d_water = ( ( (water_soil - moist_p_b) / (moist_p_a-moist_p_b) )**  &
      &         ( moist_p_d*( (moist_p_b-moist_p_a) / moist_p_a ) ) ) * &
      &       ( (water_soil / moist_p_a)**moist_p_d )

  END FUNCTION get_modifier_water

  ! ---------------------------------------------------------------------------------------------------------
  FUNCTION get_modifier_depth(layer_depth) RESULT(d_depth)
    USE mo_carbon_constants,  ONLY: ref_dep

    REAL(lp), INTENT(in)  :: layer_depth
    REAL(lp) :: d_depth

    d_depth = exp(-layer_depth/ref_dep)                       ! depth modifier as CLM

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
    REAL(lp), INTENT(out)   :: diff_coeff(:)           !< Vertical carbon mixing coefficient [m2 s-1]
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
    SUBROUTINE move_ag_2_bg(nsoil, max_depth, layer_thick, permafrost_period_sl,                          &
      &                     green_cpools_ag, woody_cpools_ag, green_cpools_bg_sl, woody_cpools_bg_sl,     &
      &                     c_humus_1_sl, c_humus_2_sl)

    USE mo_carbon_constants, ONLY: C_2_biomass_ratio, rho_org_ag, rho_org_bg, molarMassC_kg, &
      &                            permafrost_threshold

    INTEGER,  INTENT(in)    :: nsoil                   !< Number of soil layers []
    REAL(lp), INTENT(in)    :: max_depth               !< Soil depth until bedrock or permafrost [m]
    REAL(lp), INTENT(in)    :: layer_thick(:)          !< Thickness of soil layers [m]
    REAL(lp), INTENT(in)    :: permafrost_period_sl(:) !< Length of consecutive frozen ground for each soil layer [yr]
    REAL(lp), INTENT(inout) :: green_cpools_ag(:)      !< Above ground green (leaf) AWEN carbon pools [mol(C) m-2]
    REAL(lp), INTENT(inout) :: woody_cpools_ag(:)      !< Above ground woody AWEN carbon pools [mol(C) m-2]
    REAL(lp), INTENT(inout) :: green_cpools_bg_sl(:,:) !< Below ground green (leaf) AWEN carbon pools [mol(C) m-2]
    REAL(lp), INTENT(inout) :: woody_cpools_bg_sl(:,:) !< Below ground woody AWEN carbon pools [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_1_sl(:)         !< Humus pool from green litter pool [mol(C) m-2]
    REAL(lp), INTENT(inout) :: c_humus_2_sl(:)         !< Humus pool from woody litter pool [mol(C) m-2]

    REAL(lp),  PARAMETER :: thresh_frac_org    = 1.25_lp  !< Threshold organic fraction above which C is being moved downwards
    REAL(lp),  PARAMETER :: ag_layer_depth     = 0.13_lp  !< Assumed depth of above ground litter layer [m]
    REAL(lp) :: c_acid_ag1                  !< Above ground acid soluble green litter pool [mol(C) m-2]
    REAL(lp) :: c_acid_bg1_sl(nsoil)        !< Below ground acid soluble green litter pool [mol(C) m-2]
    REAL(lp) :: c_acid_ag2                  !< Above ground acid soluble woody litter pool [mol(C) m-2]
    REAL(lp) :: c_acid_bg2_sl(nsoil)        !< Below ground acid soluble woody litter pool [mol(C) m-2]
    REAL(lp) :: c_water_ag1                 !< Above ground water soluble green litter pool [mol(C) m-2]
    REAL(lp) :: c_water_bg1_sl(nsoil)       !< Below ground water soluble green litter pool [mol(C) m-2]
    REAL(lp) :: c_water_ag2                 !< Above ground water soluble woody litter pool [mol(C) m-2]
    REAL(lp) :: c_water_bg2_sl(nsoil)       !< Below ground water soluble woody litter pool [mol(C) m-2]
    REAL(lp) :: c_ethanol_ag1               !< Above ground ethanol soluble green litter pool [mol(C) m-2]
    REAL(lp) :: c_ethanol_bg1_sl(nsoil)     !< Below ground ethanol soluble green litter pool [mol(C) m-2]
    REAL(lp) :: c_ethanol_ag2               !< Above ground ethanol soluble woody litter pool [mol(C) m-2]
    REAL(lp) :: c_ethanol_bg2_sl(nsoil)     !< Below ground ethanol soluble woody litter pool [mol(C) m-2]
    REAL(lp) :: c_nonsoluble_ag1            !< Above ground nonsoluble green litter pool [mol(C) m-2]
    REAL(lp) :: c_nonsoluble_bg1_sl(nsoil)  !< Below ground nonsoluble green litter pool [mol(C) m-2]
    REAL(lp) :: c_nonsoluble_ag2            !< Above ground nonsoluble woody litter pool [mol(C) m-2]
    REAL(lp) :: c_nonsoluble_bg2_sl(nsoil)  !< Below ground nonsoluble woody litter pool [mol(C) m-2]
    REAL(lp) :: cpool_tot
    REAL(lp) :: depth_upper
    REAL(lp) :: rfrac
    REAL(lp) :: tmp

    REAL(lp) :: max_C_ag
    REAL(lp) :: max_C_bg(size(layer_thick))
    INTEGER  :: is, isp1


    ! Above ground
    c_acid_ag1             = green_cpools_ag(Y_ACID)
    c_water_ag1            = green_cpools_ag(Y_WATER)
    c_ethanol_ag1          = green_cpools_ag(Y_ETHANOL)
    c_nonsoluble_ag1       = green_cpools_ag(Y_NONSOLUBLE)
    c_acid_ag2             = woody_cpools_ag(Y_ACID)
    c_water_ag2            = woody_cpools_ag(Y_WATER)
    c_ethanol_ag2          = woody_cpools_ag(Y_ETHANOL)
    c_nonsoluble_ag2       = woody_cpools_ag(Y_NONSOLUBLE)
    ! ! Below ground
    c_acid_bg1_sl(:)       = green_cpools_bg_sl(:,Y_ACID)
    c_water_bg1_sl(:)      = green_cpools_bg_sl(:,Y_WATER)
    c_ethanol_bg1_sl(:)    = green_cpools_bg_sl(:,Y_ETHANOL)
    c_nonsoluble_bg1_sl(:) = green_cpools_bg_sl(:,Y_NONSOLUBLE)
    c_acid_bg2_sl(:)       = woody_cpools_bg_sl(:,Y_ACID)
    c_water_bg2_sl(:)      = woody_cpools_bg_sl(:,Y_WATER)
    c_ethanol_bg2_sl(:)    = woody_cpools_bg_sl(:,Y_ETHANOL)
    c_nonsoluble_bg2_sl(:) = woody_cpools_bg_sl(:,Y_NONSOLUBLE)

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

    ! Backmap to AWEN vectors
    ! Above ground
    green_cpools_ag(Y_ACID)            = c_acid_ag1
    green_cpools_ag(Y_WATER)           = c_water_ag1
    green_cpools_ag(Y_ETHANOL)         = c_ethanol_ag1
    green_cpools_ag(Y_NONSOLUBLE)      = c_nonsoluble_ag1
    woody_cpools_ag(Y_ACID)            = c_acid_ag2
    woody_cpools_ag(Y_WATER)           = c_water_ag2
    woody_cpools_ag(Y_ETHANOL)         = c_ethanol_ag2
    woody_cpools_ag(Y_NONSOLUBLE)      = c_nonsoluble_ag2
    ! ! Below ground
    green_cpools_bg_sl(:,Y_ACID)       = c_acid_bg1_sl(:)
    green_cpools_bg_sl(:,Y_WATER)      = c_water_bg1_sl(:)
    green_cpools_bg_sl(:,Y_ETHANOL)    = c_ethanol_bg1_sl(:)
    green_cpools_bg_sl(:,Y_NONSOLUBLE) = c_nonsoluble_bg1_sl(:)
    woody_cpools_bg_sl(:,Y_ACID)       = c_acid_bg2_sl(:)
    woody_cpools_bg_sl(:,Y_WATER)      = c_water_bg2_sl(:)
    woody_cpools_bg_sl(:,Y_ETHANOL)    = c_ethanol_bg2_sl(:)
    woody_cpools_bg_sl(:,Y_NONSOLUBLE) = c_nonsoluble_bg2_sl(:)

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

  ! ---------------------------------------------------------------------------------------------------------
  !>
  !> Sum up layered carbon storages for single below ground pools
  !>
  !> This routine accumulates the vertically distinct carbon fields of the multi-layer
  !> yasso scheme over all soil layers to fill the arrays used for the single
  !> layer yasso.
  !>
  SUBROUTINE sum_up_carbon_layers(nc, nsoil, cpool_sl, cpool_bg)

    INTEGER,  INTENT(in)     :: nc                    !< Number of grid cells [ ]
    INTEGER,  INTENT(in)     :: nsoil                 !< Number of soil layers [ ]
    REAL(lp), INTENT(in)     :: cpool_sl(:,:)         !< Soil carbon on layers [mol(C) m-2]
    REAL(lp), INTENT(out)    :: cpool_bg(:)        !< Below ground soil carbon pool [mol(C) m-2]

    INTEGER  :: ic                                    !< Grid cell index
    INTEGER  :: is                                    !< Soil layer index

    ! First set below ground carbon pool to zero
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      cpool_bg(ic) = 0._lp
    END DO
    !$ACC END PARALLEL LOOP

    ! Add carbon from each layer to below ground pool
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(acc_stream)
    DO ic = 1, nc
      DO is = 1, nsoil
        cpool_bg(ic) = cpool_bg(ic) * cpool_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE sum_up_carbon_layers
#endif
END MODULE mo_carbon_process
