!> Contains the routines for the carbon processes
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
MODULE mo_sb_carbon_process

  USE mo_jsb_control,      ONLY: acc_stream
  USE mo_kind,             ONLY: wp
  USE mo_lnd_bgcm_idx
  USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround, &
                               & days_per_year, sec_per_year, sec_per_day

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_jsb_C_litter, prep_jsb_sb_processes, get_per_tile

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_carbon_process'

CONTAINS

  ! Calculations to allocate NPP to the plant C pools
  ! Step 1: allocate NPP to SB litter pools
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_jsb_C_litter( &
    & LAI,                        & ! in
    & LAI_previous,               & ! in
    & max_LAI,                    & ! in
    & NPP_pot_yDayMean,           & ! in
    !
    & fract_npp_2_woodPool,       & ! in
    & fract_npp_2_reservePool,    & ! in
    & fract_npp_2_exudates,       & ! in
    & fract_green_2_herbivory,    & ! in
    !
    & tau_c_woods,                & ! in
    !
    & LAI_shed_constant,          & ! in
    & Max_C_content_woods,        & ! in
    & specific_leaf_area_C,       & ! in
    & reserveC2leafC,             & ! in
    & pheno_type,                 & ! in
    !
    & c_green,                    & ! inout
    & c_woods,                    & ! inout
    & c_reserve,                  & ! inout
    & c_crop_harvest,             & ! inout
    !
    & NPP_flux_correction,        & ! out
    & excess_NPP,                 & ! out
    & root_exudates,              & ! out
    & cflux_c_greenwood_2_litter, & ! out
    & cflux_herbivory,            & ! out
    & cflux_herb_2_littergreen,   & ! out
    & cflux_herb_2_atm,           & ! out
    & NPP_act,                    & ! out
    & leafLitter,                 & ! out
    & C_2_litter_woodPool         & ! out
    & )

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY : tmelt     ! tmelt =273.15_wp
    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    ! Use of process memory
    dsl4jsb_Use_memory(SB_)    

    ! Arguments
    REAL(wp), INTENT(in)    :: LAI                    ! Yesterdays mean LAI
    REAL(wp), INTENT(in)    :: LAI_previous           ! The day before yesterdays mean LAI
    REAL(wp), INTENT(in)    :: max_LAI                ! Maximum value of LAI
    REAL(wp), INTENT(in)    :: NPP_pot_yDayMean       ! Yesterdays mean NPP rate [mol(C)/m^2 s]
    REAL(wp), INTENT(in)    :: fract_npp_2_woodPool    ! Fraction of NPP to be put maximally into the green pool
    REAL(wp), INTENT(in)    :: fract_npp_2_reservePool ! Fraction of NPP to be put into the optimally into the reserve pool
    REAL(wp), INTENT(in)    :: fract_npp_2_exudates    ! Fraction of NPP to be put into the optimally into the root exudates
    REAL(wp), INTENT(in)    :: fract_green_2_herbivory
    REAL(wp), INTENT(in)    :: tau_c_woods        ! Time constant by which woods Pool is depreciated  !! [days] !!
                                                 ! Note, in the lct lib it is given in [years] but when it is read in
                                                 ! by JS4 it is converted to [days].
    REAL(wp), INTENT(in)    :: LAI_shed_constant      ! Leaf shedding at a constant rate for evergreens etc. [1/day]
    REAL(wp), INTENT(in)    :: Max_C_content_woods    ! Maximum carbon content of wood pool (from lctLib)
    REAL(wp), INTENT(in)    :: specific_leaf_area_C   ! Specific leaf area (from lctLib) [m^2(leaf)/mol(Carbon)]
    REAL(wp), INTENT(in)    :: reserveC2leafC         ! Ratio of max. carbon in reserve pool to max. carbon of leaves
    INTEGER,  INTENT(in)    :: pheno_type             ! crop =5
    !
    REAL(wp), INTENT(inout) :: c_green            ! Green carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_woods            ! Wood carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_reserve          ! Reserve carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_crop_harvest     ! Crop harvest carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    !REAL(wp), INTENT(out)   :: soilResp_rate_pot      ! soil respiration without N limitation
    REAL(wp), INTENT(out)   :: NPP_flux_correction    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                                     !    correction arises either because otherwise the reserve pool would get
                                                     !    negative (positive correction), or the wood pool would exceed its
                                                     !    maximum value (negative correction). [mol(C)/m^2 s]
    REAL(wp), INTENT(out)   :: excess_NPP             ! Part of NPP that could not be stored in one of the plant carbon pools
                                                     !    (green, wood, reserve), but had to be  thrown away (a posteriori
                                                     !    reduction of NPP) [mol(C)/m^2 s]
    REAL(wp), INTENT(out)   :: root_exudates          ! Value of NPP_2_rootExudates had to be dropped into the litter_green_bg pool.
                                                     !    [mol(C)/m^2 s]
    REAL(wp), INTENT(out)   :: cflux_c_greenwood_2_litter       ! Total carbon flux from vegetation to litter (green+woody, ag+bg)
    REAL(wp), INTENT(out)   :: cflux_herbivory
    REAL(wp), INTENT(out)   :: cflux_herb_2_littergreen
    REAL(wp), INTENT(out)   :: cflux_herb_2_atm
    REAL(wp), INTENT(out)   :: NPP_act               ! Actual NPP after N-limitation and excess carbon drop [mol(C)/m^2 s]
    REAL(wp), INTENT(out)   :: leafLitter            ! amount of leafLitter
    REAL(wp), INTENT(out)   :: C_2_litter_woodPool   ! amount of woodyLitter
    !
    ! added for thonicke fire algorithm (spitfire) needed to correctly split the carbon pools into the different fuel classes
    !real(wp),optional,intent(inout) :: fract_litter_wood_new   ! new fraction in above ground wood litter pool


    ! locals
    LOGICAL  :: is_crop                 ! logical crop mask

    REAL(wp) :: c_green_max,c_reserve_optimal
    REAL(wp) :: NPP_2_greenPool, NPP_2_woodPool, NPP_2_reservePool, NPP_2_rootExudates, Green_2_herbivory
    REAL(wp) :: C_2_litter_greenPool
    REAL(wp) :: excess_carbon
    REAL(wp) :: leaf_shedding

    REAL(wp) :: c_green_pot, c_woods_pot         !! values of potential carbon pools i.e. before accounting for N-limitation
    REAL(wp) :: c_reserve_pot

    REAL(wp) :: Cflx_NPP_2_green_pot,Cflx_NPP_2_wood_pot,Cflx_NPP_2_reserve_pot
    REAL(wp) :: Cflx_2_crop_harvest, Cflx_crop_harvest_2_atmos
    REAL(wp) :: Cflx_faeces_2_LG, Cflx_faeces_2_atm
    ! R: only for spitfire
    !REAL(wp) :: litter_wood_new, litter_wood_total

    ! local variables


    !REAL(wp) :: Nflux_NPP_2_green_pot,Nflux_NPP_2_wood_pot
    !REAL(wp) :: Nflux_litterGreen_ag_2_slow_pot, Nflux_litterGreen_bg_2_slow_pot
    !REAL(wp) :: Nflux_litter_wood_ag_2_slow_pot,Nflux_litter_wood_bg_2_slow_pot
    !REAL(wp) :: minNflux_litter_green_ag_pot,minNflux_litter_green_bg_pot
    !REAL(wp) :: minNflux_litter_wood_ag_pot,minNflux_litter_wood_bg_pot
    !REAL(wp) :: N2O_emissions_litter_pot

    !REAL(wp) :: minN_plant_demand,minN_soil_demand,total_minN_demand
    REAL(wp) :: NPP_2_green_via_mobileN, NPP_2_wood_via_mobileN
    !REAL(wp) :: Nmobile_2_green, Nmobile_2_wood, Nmobile_2_green_pot, Nmobile_2_wood_pot
    !REAL(wp) :: redFact_Nmobile
    REAL(wp) :: redFact_Nlimitation
    !REAL(wp) :: Nflx_faeces_2_SMINN, minNflux_litter_active
    !REAL(wp) :: Nflx_crop_harvest_2_SMINN

    !REAL(wp), PARAMETER :: N2O_ef_grazing = 0._wp  !! N2O grazing switched off; reasonable value for the emission factor: 0.0125_wp


    ! R: Normally parameters should be in mo_carbon_constants, but the following parameters are used only in this subroutine!
    !    Therefore I put them here. Is there any not structural advantage to put them into mo_carbon_constants?
    ! Parameter for Subroutine calc_jsb_C_litter
    REAL(wp), PARAMETER :: Q10                    = 1.8_wp          ! Empirical parameter for temperature dependence of
                                                                    !   heterotrophic respiration rate
    REAL(wp), PARAMETER :: kappa                  = 1.0_wp          ! Empirical parameter for soil moisture dependence of
                                                                    !   heterotrophic respiration rate
    REAL(wp), PARAMETER :: referenceTemp_Q10      = Tmelt           ! Reference temperature in the Q10 formula [Kelvin]
    REAL(wp), PARAMETER :: tau_c_reserve          = 1.0_wp *365._wp ! time constant by which reserve pool is depreciated [days]

    REAL(wp), PARAMETER :: tau_c_crop_harvest     = 1._wp*365._wp   ! time constant by which the crop harvest pool decays [days]
    REAL(wp), PARAMETER :: greenC2leafC           = 4.0_wp          ! ratio of carbon green pool (leaves, fine roots, starches,
                                                                    !    sugars) to leaf carbon
    REAL(wp), PARAMETER :: fract_C_faeces2_LG     = 0.3_wp          ! fraction of C from faces put into the litter green pool
    REAL(wp), PARAMETER :: fract_C_crop_harvest   = 0.5_wp          ! fraction of C from ag litter flux put to the crop harvest pool

    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    !soilResp_rate_pot        = 0.0_wp
    NPP_act                  = 0.0_wp
    ! redFact_Nlimitation      = 1.0_wp

    ! R: only for spitfire:
    !IF (PRESENT(fract_litter_wood_new)) THEN
    !   litter_wood_new = 0.0_wp
    !END IF

    ! Preparations
    !-------------------------------------------

    IF (pheno_type == 5) THEN
       is_crop = .TRUE.
    ELSE
       is_crop = .FALSE.
    END IF

    ! Initializations
    NPP_flux_correction        = 0.0_wp
    c_green_max                = 0.0_wp
    c_reserve_optimal          = 0.0_wp
    excess_NPP                 = 0.0_wp
    root_exudates              = 0.0_wp
    cflux_c_greenwood_2_litter = 0.0_wp
    cflux_herbivory            = 0.0_wp
    cflux_herb_2_littergreen   = 0.0_wp
    cflux_herb_2_atm           = 0.0_wp
    Cflx_2_crop_harvest        = 0.0_wp
    leafLitter                 = 0.0_wp           ! amount of leafLitter
    C_2_litter_woodPool        = 0.0_wp

    ! on all veg tiles:

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

    leaf_shedding = MAX(LAI * LAI_shed_constant,LAI_previous - LAI)  !! Leaf shedding is assured at a minimum constant rate
                                                                     !! which maybe exceeded by the decrease in LAI.

    C_2_litter_greenPool =  &                                     !! Maximally the whole carbon from the green pool is transferred
       MIN(greenC2leafC * leaf_shedding / specific_leaf_area_C, & !!    by leaf and root shedding to the green litter pool
           c_green)
    c_green = c_green - C_2_litter_greenPool              !! This removes carbon from the green pool

    excess_carbon = max(0.0_wp,c_green - c_green_max)               !! If by the reduction of LAI the maximum value of the green
                                                                    !!    pool is still smaller than the current value, there is
    IF (excess_carbon > 0._wp) THEN
      C_2_litter_greenPool = C_2_litter_greenPool + excess_carbon   !!    excess C that has also to be shedded to the green litter
      c_green = c_green - excess_carbon                             !!    pool and also subtracted from the green pool
    END IF

    IF (is_crop) THEN
      !!  Move shedded Carbon pro rata to crop harvest flux
      Cflx_2_crop_harvest = fract_green_aboveGround * fract_C_crop_harvest * C_2_litter_greenPool
    END IF

    cflux_c_greenwood_2_litter = C_2_litter_greenPool

    !! ---------------------------------------------------------------------------------
    !! Herbivory loss from c_green & Npool_green w.r.t PFTs in lctlib

    Green_2_herbivory  = c_green * fract_green_2_herbivory             !! Grazing flux indirect comput.from NPP losses
    c_green        = c_green - Green_2_herbivory                  !! Loss due to grazing

    Cflx_faeces_2_LG  = Green_2_herbivory * fract_C_faeces2_LG             !! Cflux_faeces_2_LG to be stored in litter green pool
                                                                           !!      (not considered in cflux_c_greenwood_2_litter)
                                                                           !!      (no life times (tau) enter Cflx_faeces_2_LG!)
    Cflx_faeces_2_atm = Green_2_herbivory * (1.0_wp - fract_C_faeces2_LG)  !! Cflux_faeces_2_atm to be lost to atmosphere
                                                                           !!      (dk: not added to CO2?)
    !! ---------------------------------------------------------------------------------
    !! 2. Wood shedding: Transfer of C from wood pool to wood litter pools
    !!
                                                            !! Assuming that forests continuously die at the inverse lifetime of trees
    C_2_litter_woodPool = c_woods / tau_c_woods             ! The bigger tau gets, the weaker the C pools will be depreciated:
                                                            ! tau_c_woods: [days] -> [timestep]
                                                            !! .. the shedded wood (MAX FUNCTION FOR NONWOODY PFTS)
    ! R: only for spitfire
    ! if (present(fract_litter_wood_new))     litter_wood_new = fract_wood_aboveGround * C_2_litter_woodPool
    c_woods = c_woods -  C_2_litter_woodPool                 !! .. and then all is subtracted from the wood pool.

    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter + C_2_litter_woodPool

    !! ----------------------------------------------------------------------------
    !! 3. Depletion of reserve pool: Transfer of C from reserve pool to green litter pool
    !!    Some organic carbon of the reserve pool is always lost with mortality of plants.
    !!    Note that the reserve pool contains no Nitrogen (starches and sugar are free of N).
    !!    Note also that the green litter pool has no fixed C/N-ratio (therefore no compensation
    !!    flux is needed)

    C_2_litter_greenPool = c_reserve / tau_c_reserve
    IF (is_crop) THEN
      Cflx_2_crop_harvest = Cflx_2_crop_harvest + &
          fract_green_aboveGround * fract_C_crop_harvest * C_2_litter_greenPool
    END IF
    c_reserve = c_reserve - C_2_litter_greenPool
    c_crop_harvest = c_crop_harvest + Cflx_2_crop_harvest

    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter + C_2_litter_greenPool

    !!
    !! ============== START OF POTENTIAL PLANT-CARBON ALLOCATION TO ESTIMATE PLANT-N-DEMAND =================================
    !!                       (exception: case NPP<0 handles actual allocation)

    c_green_pot   = c_green
    c_woods_pot   = c_woods
    c_reserve_pot = c_reserve

    !! Preliminary determination of distribution of NPP to the various pools (may be corrected afterwards according to
    !!    avaialable Nitrogen)

    IF (NPP_pot_yDayMean <= 0.0_wp) THEN                    !! In case of negative or zero NPP

      NPP_2_reservePool =  NPP_pot_yDayMean * sec_per_day   !! .. the plant looses C and
      c_reserve     = c_reserve + NPP_2_reservePool         !! .. we try to take it from the reserve pool
      IF (c_reserve < 0.0_wp) THEN                          !! But if thereby the reserve pool gets negative
        NPP_2_reservePool = NPP_2_reservePool - c_reserve   !! .. we absorb negative NPP only up to zero reserve pool and
        NPP_flux_correction = -c_reserve / sec_per_day      !! .. give the deficit as a flux correction back to the
        c_reserve = 0.0_wp                                  !! .. calling routine. Hence all Carbon from the reserve pool
      END IF                                                !! .. is used up.

      !! For correct handling of nitrogen limitation below (which handles only the case of positive NPP) set NPP to zero:
      NPP_2_reservePool = 0.0_wp
      NPP_2_greenPool   = 0.0_wp   !! All other transfer rates are zero
      NPP_2_woodPool    = 0.0_wp
      NPP_2_rootExudates= 0.0_wp
      c_reserve_pot = c_reserve    !! remember modified reserve pool

    ELSE  !! NPP is positive

      NPP_2_woodPool    = fract_npp_2_woodPool * NPP_pot_yDayMean * sec_per_day  !! NPP it is distributed according to predefined
                                                                                 !!    relative fraction
      NPP_2_reservePool = fract_npp_2_reservePool * NPP_pot_yDayMean * sec_per_day
      NPP_2_rootExudates= fract_npp_2_exudates * NPP_pot_yDayMean * sec_per_day
      NPP_2_greenPool   = (1.0_wp - fract_npp_2_woodPool - fract_npp_2_reservePool                 &
           - fract_npp_2_exudates) * NPP_pot_yDayMean * sec_per_day

      !! Growth of Wood Pool (structural carbon of living plants)

      c_woods_pot = c_woods_pot + NPP_2_woodPool                  !! Then it is attempted to put all NPP share into it
      excess_carbon = max(0.0_wp,c_woods_pot-Max_C_content_woods) !! .. but thereby the pool may get too large
      IF (excess_carbon > 0._wp) THEN
        c_woods_pot = c_woods_pot - excess_carbon                 !! .. so that the excess carbon has once more to be
                                                                  !! .. subtracted
        NPP_2_greenPool = NPP_2_greenPool + excess_carbon         !! .. and instead made available to the green pool
        NPP_2_woodPool = NPP_2_woodPool - excess_carbon           !! .. and the actual amount of NPP put to the wood pool
      END IF
                                                                  !! .. is less
      !! Growth of Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
      IF (c_reserve_pot < c_reserve_optimal) THEN                 !! .. If the reserve pool is smaller than optimal,
        c_reserve_pot = c_reserve_pot +  NPP_2_reservePool        !! .... it is filled by the available NPP
        excess_carbon =                                         & !! .... Thereby it may happen that it gets larger than
             max(0.0_wp,c_reserve_pot - c_reserve_optimal)        !! .... optimal so that there is excess carbon
        IF (excess_carbon > 0._wp) THEN
          c_reserve_pot = c_reserve_pot - excess_carbon           !! .... that needs not be taken up
          NPP_2_greenPool = NPP_2_greenPool + excess_carbon       !! .... but can better be used to increase the green
          NPP_2_reservePool = NPP_2_reservePool - excess_carbon   !! .... pool and the actual amount of NPP put to the
        END IF
                                                                  !! .... reserve pool is less.
      ELSE                                                        !! .... Otherwise (reserve pool is larger than optimal)
        NPP_2_greenPool = NPP_2_greenPool + NPP_2_reservePool     !! .... all NPP is left for the green pool
        NPP_2_reservePool = 0.0_wp                                !! .... so that nothing is stored in the reserve pool.
      END IF

      !! Growth of Green Pool (leaves and fine roots): (in case of too much NPP, try to put it into the reserve pool).

      c_green_pot = c_green_pot + NPP_2_greenPool                 !! Green pool is filled by the available NPP.
      excess_carbon = max(0.0_wp,c_green_pot - c_green_max)       !! .. Thereby it may get larger than appropriate for
                                                                  !!    current LAI.
      IF (excess_carbon > 0._wp) THEN
        NPP_2_greenPool = NPP_2_greenPool - excess_carbon         !! .. Hence the actual amount of NPP put to the green
                                                                  !!    pool is less.
        c_green_pot = c_green_pot - excess_carbon                 !! .. This excess carbon needs not be taken up but
      END IF
      IF (c_reserve_pot < c_reserve_optimal) THEN                 !! .... if the reserve pool is smaller than optimal
        IF (excess_carbon > 0._wp) THEN
          c_reserve_pot = c_reserve_pot + excess_carbon           !! .... it is tried to put the carbon there,
          NPP_2_reservePool = NPP_2_reservePool + excess_carbon   !! .... which means that additional NPP is put to the
                                                                  !!      reserve pool.
        END IF
        excess_carbon =                                        &  !! .... Thereby it may happen that the reserve pool
             max(0.0_wp,c_reserve_pot - c_reserve_optimal)        !!      increases beyond the optimal value.
        IF (excess_carbon > 0._wp) THEN
          c_reserve_pot = c_reserve_pot - excess_carbon           !! .... In that case the excess carbon is once more
                                                                  !!      removed from the reserve pool
          NPP_2_reservePool = NPP_2_reservePool - excess_carbon   !! .... so that the actual amount of NPP put into the
                                                                  !!      reserve pool is less.
        END IF
      END IF

      excess_NPP = excess_carbon / sec_per_day

      root_exudates = NPP_2_rootExudates / sec_per_day

    END IF !! NPP > 0 end

    !! ============== START OF POTENTIAL LITTER-CARBON ALLOCATION ==============================================

    !! ============== use nitrogen from N mobile pool to satisfy (at least partly) for potential carbon fluxes ===
    ! Remarks: Nmobile_2_green + Nmobile_2_wood <= Npool_mobile
    ! Npool_mobile accounts minN_plant_demand by trasfering N into green and wood pool

    Cflx_NPP_2_green_pot    = NPP_2_greenPool     !! potential growth of green carbon
    Cflx_NPP_2_wood_pot     = NPP_2_woodPool      !! potential growth of wood carbon

    NPP_2_green_via_mobileN = 0.0_wp
    NPP_2_wood_via_mobileN  = 0.0_wp

    Cflx_NPP_2_green_pot = Cflx_NPP_2_green_pot - NPP_2_green_via_mobileN
    Cflx_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot  - NPP_2_wood_via_mobileN

    ! R: These lines were calculated also with_Nitrogen=.FALSE. in JS3. However, it was not used further when with_Nitrogen=.FALSE.
    !    In JS3 it should stand within the brackets.
    !Nflux_NPP_2_green_pot = Cflx_NPP_2_green_pot/cn_green
    !Nflux_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot/cn_woods
    !minN_plant_demand = Nflux_NPP_2_green_pot + Nflux_NPP_2_wood_pot

    !! ============== COLLECT POTENTIAL C-FLUXES ===============================================================

    Cflx_NPP_2_reserve_pot        = NPP_2_reservePool                                     !! potential growth of reserve carbon

    redFact_Nlimitation = 1.0_wp   !! This assures absence of nitrogen limitation in this mixed C-N code

    !! ============== UPDATE ALL C-POOLS =======================================================================

    c_green       = c_green       + redFact_Nlimitation * Cflx_NPP_2_green_pot
    c_woods       = c_woods       + redFact_Nlimitation * Cflx_NPP_2_wood_pot
    c_reserve     = c_reserve     + Cflx_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N

    !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

    IF (NPP_pot_yDayMean > 0.0_wp) THEN
      NPP_act = (NPP_2_green_via_mobileN + NPP_2_wood_via_mobileN + NPP_2_rootExudates +               &
        & redFact_Nlimitation * (Cflx_NPP_2_green_pot + Cflx_NPP_2_wood_pot) + Cflx_NPP_2_reserve_pot  &
        & ) / sec_per_day
    ELSE
      NPP_act = NPP_pot_yDayMean + NPP_flux_correction  ! NPP_pot_yDayMean is negative: NPP_flux_correction is zero
                                                                 ! except c_reserve was < 0.0_wp.
    END IF


    !! ============== COMPUTE SOIL RESPIRATION ===================================================================
    leafLitter    = cflux_c_greenwood_2_litter - C_2_litter_woodPool + Cflx_faeces_2_LG - Cflx_2_crop_harvest

    ! always:
    !! ============== Decay of crop harvest =============

    Cflx_crop_harvest_2_atmos = c_crop_harvest / tau_c_crop_harvest
    c_crop_harvest = c_crop_harvest - Cflx_crop_harvest_2_atmos

    !! ============== COMPUTE Herbivory RESPIRATION =====

    cflux_herb_2_atm = -(Cflx_faeces_2_atm) /sec_per_day

    !! ============== Litter and herbivory diagnostics (flux from day-1 to s-1)
    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter /sec_per_day
    cflux_herbivory    = Green_2_herbivory / sec_per_day
    cflux_herb_2_littergreen = Cflx_faeces_2_LG / sec_per_day   ! DSG: goes into yasso  

  END SUBROUTINE calc_jsb_C_litter




  ! --- calc_jsb_to_sb_litter_partitioning() ---------------------------------------------------------------------------------------------------
  !
  ! partition jsb leafLitter and C_2_litter_woodPool to SB litter
  !
  SUBROUTINE prep_jsb_sb_processes( &
    & nc, &                            ! in
    & nsoil_w, &
    & num_sl_above_bedrock, &
    & lctlib_sla, &
    & lctlib_growthform, &
    & sb_model_scheme, &
    & soil_depth_sl, &
    & root_fraction_sl, &
    & veg_litterfall_mt, &             ! in
    & sb_formation_mt)
    
    
    USE mo_veg_constants,           ONLY: ITREE, IGRASS, ICROP
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_sb_constants,            ONLY: fc_soluable_max, k_fc_soluable, k_fn_soluable, k_fp_soluable, &
                                          lc_leaf_max, lc_leaf2sla, lc_fine_root, lc_coarse_root, lc_fruit, lc_seed_bed
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                        !< dimension
    INTEGER,                           INTENT(in)    :: nsoil_w                   !< nr soil layer
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock      !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp),                          INTENT(in)    :: lctlib_sla                !< lctlib parameter
    INTEGER,                           INTENT(in)    :: lctlib_growthform         !< lctlib parameter
    CHARACTER(len=*),                  INTENT(in)    :: sb_model_scheme           !< config option: soil model
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: soil_depth_sl             !< soil layer thickness
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: root_fraction_sl          !< root fraction per soil layer
    REAL(wp),                          INTENT(in)    :: veg_litterfall_mt(:,:,:)  !< bgcm flux: vegetation litterfall
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)  !< bgcm flux: soil biogeochemistry formation flux
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc) :: fc_soluable_leaf           !< fc: ...
    REAL(wp), DIMENSION(nc) :: fc_soluable_fine_root
    REAL(wp), DIMENSION(nc) :: fc_soluable_coarse_root
    REAL(wp), DIMENSION(nc) :: fc_soluable_fruit
    REAL(wp), DIMENSION(nc) :: fc_soluable_seed_bed
    REAL(wp), DIMENSION(nc) :: fc_soluable_labile
    REAL(wp), DIMENSION(nc) :: fc_soluable_reserve
    REAL(wp), DIMENSION(nc) :: fn_soluable_leaf           !< fn: ...
    REAL(wp), DIMENSION(nc) :: fn_soluable_fine_root
    REAL(wp), DIMENSION(nc) :: fn_soluable_coarse_root
    REAL(wp), DIMENSION(nc) :: fn_soluable_fruit
    REAL(wp), DIMENSION(nc) :: fn_soluable_seed_bed
    REAL(wp), DIMENSION(nc) :: fn_soluable_labile
    REAL(wp), DIMENSION(nc) :: fn_soluable_reserve
    REAL(wp), DIMENSION(nc) :: fp_soluable_leaf           !< fp: ...
    REAL(wp), DIMENSION(nc) :: fp_soluable_fine_root
    REAL(wp), DIMENSION(nc) :: fp_soluable_coarse_root
    REAL(wp), DIMENSION(nc) :: fp_soluable_fruit
    REAL(wp), DIMENSION(nc) :: fp_soluable_seed_bed
    REAL(wp), DIMENSION(nc) :: fp_soluable_labile
    REAL(wp), DIMENSION(nc) :: fp_soluable_reserve
    INTEGER                 :: ic, is                     !< loop over dimensions
    INTEGER                 :: isoil                      !< loop over soil layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_partitioning'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 init local variables
    !>  labile and reserve compartments are init with 1.0
    !>
    fc_soluable_leaf        = 0.0_wp
    fc_soluable_fine_root   = 0.0_wp
    fc_soluable_coarse_root = 0.0_wp
    fc_soluable_fruit       = 0.0_wp
    fc_soluable_seed_bed    = 0.0_wp
    fn_soluable_leaf        = 0.0_wp
    fn_soluable_fine_root   = 0.0_wp
    fn_soluable_coarse_root = 0.0_wp
    fn_soluable_fruit       = 0.0_wp
    fn_soluable_seed_bed    = 0.0_wp
    fp_soluable_leaf        = 0.0_wp
    fp_soluable_fine_root   = 0.0_wp
    fp_soluable_coarse_root = 0.0_wp
    fp_soluable_fruit       = 0.0_wp
    fp_soluable_seed_bed    = 0.0_wp
    ! labile & reserve
    fc_soluable_labile      = 1.0_wp
    fc_soluable_reserve     = 1.0_wp
    fn_soluable_labile      = 1.0_wp
    fn_soluable_reserve     = 1.0_wp
    fp_soluable_labile      = 1.0_wp
    fp_soluable_reserve     = 1.0_wp

    !>1.0 determine allocation of fresh litter fall to litter pools
    !>
    DO ic = 1,nc
      ! leaf
      IF(veg_litterfall_mt(ix_leaf, ixN, ic) > eps8) THEN
        fc_soluable_leaf(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * (lc_leaf_max + lc_leaf2sla * lctlib_sla) &
          &                    * veg_litterfall_mt(ix_leaf, ixC, ic) / veg_litterfall_mt(ix_leaf, ixN, ic))
        IF(fc_soluable_leaf(ic) > eps8) THEN
          fn_soluable_leaf(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(ic)) / (k_fn_soluable * fc_soluable_leaf(ic)))
          fp_soluable_leaf(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(ic)) / (k_fp_soluable * fc_soluable_leaf(ic)))
        END IF
      END IF
      ! fine_root
      IF(veg_litterfall_mt(ix_fine_root, ixN, ic) > eps8) THEN
        fc_soluable_fine_root(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_fine_root &
          &                         * veg_litterfall_mt(ix_fine_root, ixC, ic) / veg_litterfall_mt(ix_fine_root, ixN, ic))
        IF(fc_soluable_fine_root(ic) > eps8) THEN
          fn_soluable_fine_root(ic) = 1._wp &
            &                         / (1._wp + (1._wp - fc_soluable_fine_root(ic)) / (k_fn_soluable * fc_soluable_fine_root(ic)))
          fp_soluable_fine_root(ic) = 1._wp &
            &                         / (1._wp + (1._wp - fc_soluable_fine_root(ic)) / (k_fp_soluable * fc_soluable_fine_root(ic)))
        END IF
      END IF
      ! coarse_root
      IF(veg_litterfall_mt(ix_coarse_root, ixN, ic) > eps8) THEN
        fc_soluable_coarse_root(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_coarse_root &
          &                           * veg_litterfall_mt(ix_coarse_root, ixC, ic) / veg_litterfall_mt(ix_coarse_root, ixN, ic))
        IF(fc_soluable_coarse_root(ic) > eps8) THEN
          fn_soluable_coarse_root(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(ic)) &
            &                           / (k_fn_soluable*fc_soluable_coarse_root(ic)))
          fp_soluable_coarse_root(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(ic)) &
            &                           / (k_fp_soluable*fc_soluable_coarse_root(ic)))
         END IF
      END IF
      ! fruit
      IF(veg_litterfall_mt(ix_fruit, ixN, ic) > eps8) THEN
        fc_soluable_fruit(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_fruit &
          &                     * veg_litterfall_mt(ix_fruit, ixC, ic) / veg_litterfall_mt(ix_fruit, ixN, ic))
        IF(fc_soluable_fruit(ic) > eps8) THEN
          fn_soluable_fruit(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(ic)) / (k_fn_soluable * fc_soluable_fruit(ic)))
          fp_soluable_fruit(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(ic)) / (k_fp_soluable * fc_soluable_fruit(ic)))
        END IF
      END IF
      ! seed_bed
      IF(veg_litterfall_mt(ix_seed_bed, ixN, ic) > eps8) THEN
        fc_soluable_seed_bed(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_seed_bed &
          &                        * veg_litterfall_mt(ix_seed_bed, ixC, ic) / veg_litterfall_mt(ix_seed_bed, ixN, ic))
        IF(fc_soluable_seed_bed(ic) > eps8) THEN
          fn_soluable_seed_bed(ic) = 1._wp &
            &                        / (1._wp + (1._wp - fc_soluable_seed_bed(ic)) / (k_fn_soluable * fc_soluable_seed_bed(ic)))
          fp_soluable_seed_bed(ic) = 1._wp &
            &                        / (1._wp + (1._wp - fc_soluable_seed_bed(ic)) / (k_fp_soluable * fc_soluable_seed_bed(ic)))
        END IF
      END IF
    END DO ! loop over nc

    !>2.0 Update surface metabolic and soluable litter pools with litter fall
    !>
    !>  2.1 soluable litter
    !>
    sb_formation_mt(ix_soluable_litter, ixC, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC, :)   &
      &                                 + fc_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixC, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixN, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixN, :, 1) &
      &                              + (fn_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixN, :)      &
      &                                 + fn_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixN, :)     &
      &                                 + fn_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixN, :)    &
      &                                 + fn_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixN, :)   &
      &                                 + fn_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixN, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixP, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixP, :, 1) &
      &                              + (fp_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixP, :)      &
      &                                 + fp_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixP, :)     &
      &                                 + fp_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixP, :)    &
      &                                 + fp_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixP, :)   &
      &                                 + fp_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixP, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixC13, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC13, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC13, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC13, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC13, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC13, :)   &
      &                                 + fc_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixC13, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixC14, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC14, :, 1) &
      &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC14, :)      &
      &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC14, :)     &
      &                                 + fc_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixC14, :)    &
      &                                 + fc_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixC14, :)   &
      &                                 + fc_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixC14, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixN15, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixN15, :, 1) &
      &                              + (fn_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixN15, :)      &
      &                                 + fn_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixN15, :)     &
      &                                 + fn_soluable_labile(:)   * veg_litterfall_mt(ix_labile, ixN15, :)    &
      &                                 + fn_soluable_reserve(:)  * veg_litterfall_mt(ix_reserve, ixN15, :)   &
      &                                 + fn_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixN15, :)) &
      &                              / soil_depth_sl(:,1)
    !>  2.2 polymeric litter
    !>
    sb_formation_mt(ix_polymeric_litter, ixC, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC, :, 1) &
      &                              + ((1.0_wp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC, :)      &
      &                                 + (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC, :)     &
      &                                 + (1.0_wp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC, :)    &
      &                                 + (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC, :)   &
      &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixC, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixN, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixN, :, 1) &
      &                              + ((1.0_wp - fn_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixN, :)      &
      &                                 + (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixN, :)     &
      &                                 + (1.0_wp - fn_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixN, :)    &
      &                                 + (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixN, :)   &
      &                                 + (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixN, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixP, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixP, :, 1) &
      &                              + ((1.0_wp - fp_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixP, :)      &
      &                                 + (1.0_wp - fp_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixP, :)     &
      &                                 + (1.0_wp - fp_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixP, :)    &
      &                                 + (1.0_wp - fp_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixP, :)   &
      &                                 + (1.0_wp - fp_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixP, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) &
      &                              + ((1.0_wp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC13, :)      &
      &                                 + (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC13, :)     &
      &                                 + (1.0_wp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC13, :)    &
      &                                 + (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC13, :)   &
      &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixC13, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) &
      &                              + ((1.0_wp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC14, :)      &
      &                                 + (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC14, :)     &
      &                                 + (1.0_wp - fc_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixC14, :)    &
      &                                 + (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixC14, :)   &
      &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixC14, :)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) &
      &                              + ((1.0_wp - fn_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixN15, :)      &
      &                                 + (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixN15, :)     &
      &                                 + (1.0_wp - fn_soluable_labile(:))   * veg_litterfall_mt(ix_labile, ixN15, :)    &
      &                                 + (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall_mt(ix_reserve, ixN15, :)   &
      &                                 + (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixN15, :)) &
      &                              / soil_depth_sl(:,1)

    !>3.0 Add surface litter fall from sapwood and heartwood pools to surface woody litter pool for trees
    !>    for grasses, stems and halms (=sapwood) are not actually wood, and therefore paritioned to
    !>    metabolic and structural litter as if their were fine roots
    SELECT CASE(lctlib_growthform)
    ! tree PFT
    CASE (ITREE)
      ! woody litter
      ! @TODO may be re-written to a one-liner using bgcm matrices
      sb_formation_mt(ix_woody_litter, ixC, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixN, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixN, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixP, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixP, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixC13, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC13, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixC14, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC14, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixN15, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixN15, :, 1) &
        &           + (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) / soil_depth_sl(:,1)
    ! grass PFT
    CASE (IGRASS,ICROP)
      ! soluable litter
      sb_formation_mt(ix_soluable_litter, ixC, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixN, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixN, :, 1) &
        &             + fn_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixP, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixP, :, 1) &
        &             + fp_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixC13, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC13, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixC14, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC14, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixN15, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixN15, :, 1) &
        &             + fn_soluable_fine_root(:) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) &
        &             / soil_depth_sl(:,1)
      ! polymeric litter
      sb_formation_mt(ix_polymeric_litter, ixC, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC, :, 1) &
        &             + (1.0_wp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC, :) + veg_litterfall_mt(ix_heart_wood, ixC, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixN, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixN, :, 1) &
        &             + (1.0_wp - fn_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN, :) + veg_litterfall_mt(ix_heart_wood, ixN, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixP, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixP, :, 1) &
        &             + (1.0_wp - fp_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixP, :) + veg_litterfall_mt(ix_heart_wood, ixP, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) &
        &             + (1.0_wp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) &
        &             + (1.0_wp - fc_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) &
        &             + (1.0_wp - fn_soluable_fine_root(:)) &
        &             * (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) &
        &             / soil_depth_sl(:,1)
    END SELECT

    !>4.0 Distribute fine and coarse root litter to soil depth according to the root profile
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        ! soluable litter (fine root & coarse root)
        sb_formation_mt(ix_soluable_litter, ixC, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixN, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixN, ic, is) &
          &             + (fn_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixN, ic) &
          &             +  fn_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixP, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixP, ic, is) &
          &             + (fp_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixP, ic) &
          &             +  fp_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixC13, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC13, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC13, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixC14, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC14, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC14, ic) &
          &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixN15, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixN15, ic, is) &
          &             + (fn_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixN15, ic) &
          &             +  fn_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        ! simple soil model / jsm
        SELECT CASE(TRIM(sb_model_scheme))
        CASE ("simple_1d")
          ! polymeric litter (fine root & coarse root)
          sb_formation_mt(ix_polymeric_litter, ixC, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC, ic) &
            &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixN, ic) &
            &           +  (1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixP, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixP, ic) &
            &           +  (1.0_wp - fp_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC13, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC13, ic) &
            &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC14, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC14, ic) &
            &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN15, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixN15, ic) &
            &           +  (1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
        CASE ("jsm")
          ! polymeric litter (fine root)
          sb_formation_mt(ix_polymeric_litter, ixC, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixP, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! polymeric litter (coarse root)
          sb_formation_mt(ix_woody_litter, ixC, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixN, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixP, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixP, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixC13, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC13, ic, is) &
            &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixC14, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC14, ic, is) &
            &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixN15, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixN15, ic, is) &
            &           + ((1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
        END SELECT
      END DO
    END DO

  END SUBROUTINE prep_jsb_sb_processes
  
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
  
END MODULE mo_sb_carbon_process
