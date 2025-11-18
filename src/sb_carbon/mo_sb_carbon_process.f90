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

  PUBLIC :: calc_jsb_C_litter, prep_jsb_sb_processes, get_per_tile, convert_lct_to_sla

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
    & c_leaf,                     & ! inout
    & c_fineroot,                 & ! inout
    & c_fruit,                    & ! inout
    & c_seedbed,                  & ! inout
    & c_sapheartwood,             & ! inout
    & c_coarseroot,               & ! inout
    !
    & NPP_flux_correction,        & ! out
    & excess_NPP,                 & ! out
    & root_exudates,              & ! out
    & cflux_c_greenwood_2_litter, & ! out
    & cflux_herbivory,            & ! out
    & cflux_herb_2_littergreen,   & ! out
    & cflux_herb_2_atm,           & ! out
    & NPP_act,                    & ! out
    & litterfall_leaf,            & ! out
    & litterfall_fineroot,        & ! out
    & litterfall_fruit,           & ! out
    & litterfall_seedbed,         & ! out
    & litterfall_wood,            & ! out
    & litterfall_coarseroot       & ! out
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
    REAL(wp), INTENT(inout) :: c_leaf             ! Green carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_fineroot         ! Wood carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_fruit            ! Reserve carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_seedbed          ! Crop harvest carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_sapheartwood     ! Reserve carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_coarseroot       ! Crop harvest carbon pool: on input last value; updated on output
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

    REAL(wp), INTENT(inout)   :: litterfall_leaf            ! amount of leafLitter
    REAL(wp), INTENT(inout)   :: litterfall_fineroot   ! amount of woodyLitter
    REAL(wp), INTENT(inout)   :: litterfall_fruit            ! amount of leafLitter
    REAL(wp), INTENT(inout)   :: litterfall_seedbed   ! amount of woodyLitter
    REAL(wp), INTENT(inout)   :: litterfall_wood            ! amount of leafLitter
    REAL(wp), INTENT(inout)   :: litterfall_coarseroot   ! amount of woodyLitter

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

    redFact_Nlimitation  = 1.0_wp   !! This assures absence of nitrogen limitation in this mixed C-N code
    Cflx_NPP_2_green_pot = redFact_Nlimitation * Cflx_NPP_2_green_pot
    Cflx_NPP_2_wood_pot  = redFact_Nlimitation * Cflx_NPP_2_wood_pot
    !! ============== COMPUTE SOIL RESPIRATION ===================================================================
    leafLitter    = cflux_c_greenwood_2_litter - C_2_litter_woodPool + Cflx_faeces_2_LG - Cflx_2_crop_harvest
    
    !! ============== PARTITION JSB C & LITTER FLUX TO QUINCY ===================================================
    ! partition leaf and woody litter to QUINCY litter types
    ! leafLitter --> leaf, fine_root, fruit, seedbed
    ! C_2_litter_woodPool --> coarse_root, and wood
    CALL quincy_allometry_turnover(pheno_type, &
         &                    leafLitter, C_2_litter_woodPool, &
         &                    c_green, c_woods, c_reserve, &
         &                    Cflx_NPP_2_green_pot, Cflx_NPP_2_wood_pot, Cflx_NPP_2_reserve_pot, &
         &                    c_leaf, c_fineroot, c_fruit, c_seedbed, c_sapheartwood, c_coarseroot, & 
         &                    litterfall_leaf, litterfall_fineroot, &
         &                    litterfall_fruit, litterfall_seedbed, &
         &                    litterfall_wood, litterfall_coarseroot)

    !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

    IF (NPP_pot_yDayMean > 0.0_wp) THEN
      NPP_act = (NPP_2_green_via_mobileN + NPP_2_wood_via_mobileN + NPP_2_rootExudates +               &
        & redFact_Nlimitation * (Cflx_NPP_2_green_pot + Cflx_NPP_2_wood_pot) + Cflx_NPP_2_reserve_pot  &
        & ) / sec_per_day
    ELSE
      NPP_act = NPP_pot_yDayMean + NPP_flux_correction  ! NPP_pot_yDayMean is negative: NPP_flux_correction is zero
                                                                 ! except c_reserve was < 0.0_wp.
    END IF




    ! always:
    !! ============== Decay of crop harvest =============

    Cflx_crop_harvest_2_atmos = c_crop_harvest / tau_c_crop_harvest
    c_crop_harvest = c_crop_harvest - Cflx_crop_harvest_2_atmos

    !! ============== COMPUTE Herbivory RESPIRATION =====

    cflux_herb_2_atm = -(Cflx_faeces_2_atm) /sec_per_day

    !! ============== Litter and herbivory diagnostics (flux from day-1 to s-1)
    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter /sec_per_day
    cflux_herbivory    = Green_2_herbivory / sec_per_day
    cflux_herb_2_littergreen = Cflx_faeces_2_LG / sec_per_day
    
    

    litterfall_leaf       = litterfall_leaf / sec_per_day
    litterfall_fineroot   = litterfall_fineroot / sec_per_day
    litterfall_fruit      = litterfall_fruit / sec_per_day
    litterfall_seedbed    = litterfall_seedbed / sec_per_day
    litterfall_wood       = litterfall_wood / sec_per_day
    litterfall_coarseroot = litterfall_coarseroot / sec_per_day


  END SUBROUTINE calc_jsb_C_litter




  ! --- calc_jsb_to_sb_litter_partitioning() ---------------------------------------------------------------------------------------------------
  !
  ! partition jsb leafLitter and C_2_litter_woodPool to SB litter
  !
  SUBROUTINE prep_jsb_sb_processes( &
    & nc, &                            ! in
    & nsoil_w, &
    & dtime, &
    & num_sl_above_bedrock, &
    & lctlib_sla, &
    & lctlib_cn_leaf, &
    & lctlib_np_leaf, &
    & lctlib_PhenologyType, &
    & sb_model_scheme, &
    & soil_depth_sl, &
    & root_fraction_sl, &
    & leafLitter, woodLitter, &             ! in
    & sb_formation_mt)
    
    
    USE mo_veg_constants,           ONLY: ITREE, IGRASS, ICROP
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_sb_constants,            ONLY: fc_soluable_max, k_fc_soluable, k_fn_soluable, k_fp_soluable, &
                                          lc_leaf_max, lc_leaf2sla, lc_fine_root, lc_coarse_root, lc_fruit, lc_seed_bed, &
                                          root2leaf_cn, wood2leaf_cn, root2leaf_np, wood2leaf_np, resorp_fract_leaf, resorp_fract_wood

    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                        !< dimension
    INTEGER,                           INTENT(in)    :: nsoil_w                   !< nr soil layer
    REAL(wp),                          INTENT(in)    :: dtime                     !< 
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock      !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp),                          INTENT(in)    :: lctlib_sla                !< 
    REAL(wp),                          INTENT(in)    :: lctlib_cn_leaf            !< 
    REAL(wp),                          INTENT(in)    :: lctlib_np_leaf            !< 
    INTEGER,                           INTENT(in)    :: lctlib_PhenologyType      !< JSB phenotype
    CHARACTER(len=*),                  INTENT(in)    :: sb_model_scheme           !< config option: soil model
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: soil_depth_sl             !< soil layer thickness
    REAL(wp), DIMENSION(nc, nsoil_w),  INTENT(in)    :: root_fraction_sl          !< root fraction per soil layer
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: leafLitter                !< jsb flux: vegetation leaf litterfall [mol C/m2 /s]
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: woodLitter                !< jsb flux: vegetation wood litterfall [mol C/m2 /s]   
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)  !< bgcm flux: soil biogeochemistry formation flux [mol C/m3/ dt]
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
    REAL(wp), DIMENSION(nc) :: litterfall_leaf, litterfall_fruit, litterfall_seedbed, litterfall_fineroot, &
                               litterfall_wood, litterfall_coarseroot  
    REAL(wp),               :: leaf_litter_cn, fineroot_litter_cn, coarseroot_litter_cn, &
                               leaf_litter_np, fineroot_litter_np, coarseroot_litter_np
    INTEGER                 :: ic, is                     !< loop over dimensions
    INTEGER                 :: isoil                      !< loop over soil layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_partitioning'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.1 init local variables
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
    !> litter cn ratio
    !> leaf and wood recycle nutrient before litterfall
    !> fruit and seed bed have the same cn as fresh leaf
    leaf_litter_cn          = lctlib_cn_leaf / resorp_fract_leaf
    fineroot_litter_cn      = lctlib_cn_leaf / root2leaf_cn
    coarseroot_litter_cn    = lctlib_cn_leaf / wood2leaf_cn
    wood_litter_cn          = lctlib_cn_leaf / wood2leaf_cn / resorp_fract_wood
    
    leaf_litter_np          = lctlib_np_leaf / resorp_fract_leaf
    fineroot_litter_np      = lctlib_np_leaf / root2leaf_np
    coarseroot_litter_np    = lctlib_np_leaf / wood2leaf_np
    wood_litter_np          = lctlib_np_leaf / wood2leaf_np / resorp_fract_wood

    !>1.0 convert unit of jsb leaf and woody litter [mol/s] to SB flux [mol/dtime]
    ! 
    litterfall_leaf = leafLitter * dtime
    litterfall_wood = woodLitter * dtime


    !>2.0 determine allocation of fresh litter fall to litter pools
    !>
    DO ic = 1,nc
      ! leaf
      fc_soluable_leaf(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * (lc_leaf_max + lc_leaf2sla * lctlib_sla) &
        &                    * leaf_litter_cn)
      IF(fc_soluable_leaf(ic) > eps8) THEN
        fn_soluable_leaf(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(ic)) / (k_fn_soluable * fc_soluable_leaf(ic)))
        fp_soluable_leaf(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(ic)) / (k_fp_soluable * fc_soluable_leaf(ic)))
      END IF
      ! fine_root
      fc_soluable_fine_root(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_fine_root &
        &                         * fineroot_litter_cn)
      IF(fc_soluable_fine_root(ic) > eps8) THEN
        fn_soluable_fine_root(ic) = 1._wp &
          &                         / (1._wp + (1._wp - fc_soluable_fine_root(ic)) / (k_fn_soluable * fc_soluable_fine_root(ic)))
        fp_soluable_fine_root(ic) = 1._wp &
          &                         / (1._wp + (1._wp - fc_soluable_fine_root(ic)) / (k_fp_soluable * fc_soluable_fine_root(ic)))
      END IF
      ! coarse_root
      fc_soluable_coarse_root(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_coarse_root &
        &                           * coarseroot_litter_cn))
      IF(fc_soluable_coarse_root(ic) > eps8) THEN
        fn_soluable_coarse_root(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(ic)) &
          &                           / (k_fn_soluable*fc_soluable_coarse_root(ic)))
        fp_soluable_coarse_root(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(ic)) &
          &                           / (k_fp_soluable*fc_soluable_coarse_root(ic)))
       END IF
      ! fruit, same cn as leaf litter
      fc_soluable_fruit(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_fruit &
        &                     * leaf_litter_cn)
      IF(fc_soluable_fruit(ic) > eps8) THEN
        fn_soluable_fruit(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(ic)) / (k_fn_soluable * fc_soluable_fruit(ic)))
        fp_soluable_fruit(ic) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(ic)) / (k_fp_soluable * fc_soluable_fruit(ic)))
      END IF
      ! seed_bed, same cn as fruit
      fc_soluable_seed_bed(ic) = MAX(0.0_wp, fc_soluable_max - k_fc_soluable * lc_seed_bed &
        &                        * leaf_litter_cn)
      IF(fc_soluable_seed_bed(ic) > eps8) THEN
        fn_soluable_seed_bed(ic) = 1._wp &
          &                        / (1._wp + (1._wp - fc_soluable_seed_bed(ic)) / (k_fn_soluable * fc_soluable_seed_bed(ic)))
        fp_soluable_seed_bed(ic) = 1._wp &
          &                        / (1._wp + (1._wp - fc_soluable_seed_bed(ic)) / (k_fp_soluable * fc_soluable_seed_bed(ic)))
      END IF
    END DO ! loop over nc

    !>3.0 Update surface metabolic and soluable litter pools with litter fall
    !>
    !>  3.1 soluable litter
    !>
    sb_formation_mt(ix_soluable_litter, ixC, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixC, :, 1) &
      &                              + (fc_soluable_leaf(:)       * litterfall_leaf      &
      &                                 + fc_soluable_fruit(:)    * litterfall_fruit     &
      &                                 + fc_soluable_seed_bed(:) * litterfall_seedbed) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixN, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixN, :, 1) &
      &                              + (fn_soluable_leaf(:)       * (litterfall_leaf / leaf_litter_cn)      &
      &                                 + fn_soluable_fruit(:)    * (litterfall_fruit / lctlib_cn_leaf)     &
      &                                 + fn_soluable_seed_bed(:) * (litterfall_seedbed / lctlib_cn_leaf)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_soluable_litter, ixP, :, 1) = &
      &                              sb_formation_mt(ix_soluable_litter, ixP, :, 1) &
      &                              + (fp_soluable_leaf(:)       * (litterfall_leaf / leaf_litter_cn / leaf_litter_np) &
      &                                 + fp_soluable_fruit(:)    * (litterfall_fruit / lctlib_cn_leaf / lctlib_np_leaf)     &
      &                                 + fp_soluable_seed_bed(:) * (litterfall_seedbed / lctlib_cn_leaf / lctlib_np_leaf)) &
      &                              / soil_depth_sl(:,1)
    ! sb_formation_mt(ix_soluable_litter, ixC13, :, 1) = &
      ! &                              sb_formation_mt(ix_soluable_litter, ixC13, :, 1) &
      ! &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC13, :)      &
      ! &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC13, :)     &
      ! &                                 + fc_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixC13, :)) &
      ! &                              / soil_depth_sl(:,1)
    ! sb_formation_mt(ix_soluable_litter, ixC14, :, 1) = &
      ! &                              sb_formation_mt(ix_soluable_litter, ixC14, :, 1) &
      ! &                              + (fc_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixC14, :)      &
      ! &                                 + fc_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixC14, :)     &
      ! &                                 + fc_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixC14, :)) &
      ! &                              / soil_depth_sl(:,1)
    ! sb_formation_mt(ix_soluable_litter, ixN15, :, 1) = &
      ! &                              sb_formation_mt(ix_soluable_litter, ixN15, :, 1) &
      ! &                              + (fn_soluable_leaf(:)       * veg_litterfall_mt(ix_leaf, ixN15, :)      &
      ! &                                 + fn_soluable_fruit(:)    * veg_litterfall_mt(ix_fruit, ixN15, :)     &
      ! &                                 + fn_soluable_seed_bed(:) * veg_litterfall_mt(ix_seed_bed, ixN15, :)) &
      ! &                              / soil_depth_sl(:,1)
    !>  3.2 polymeric litter
    !>
    sb_formation_mt(ix_polymeric_litter, ixC, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixC, :, 1) &
      &                              + ((1.0_wp - fc_soluable_leaf(:))       * litterfall_leaf      &
      &                                 + (1.0_wp - fc_soluable_fruit(:))    * litterfall_fruit     &
      &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * litterfall_seedbed) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixN, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixN, :, 1) &
      &                              + ((1.0_wp - fn_soluable_leaf(:))       * (litterfall_leaf / leaf_litter_cn)      &
      &                                 + (1.0_wp - fn_soluable_fruit(:))    * (litterfall_fruit / lctlib_cn_leaf)     &
      &                                 + (1.0_wp - fn_soluable_seed_bed(:)) * (litterfall_seedbed / lctlib_cn_leaf)) &
      &                              / soil_depth_sl(:,1)
    sb_formation_mt(ix_polymeric_litter, ixP, :, 1) = &
      &                              sb_formation_mt(ix_polymeric_litter, ixP, :, 1) &
      &                              + ((1.0_wp - fp_soluable_leaf(:))       * (litterfall_leaf / leaf_litter_cn / leaf_litter_np) &
      &                                 + (1.0_wp - fp_soluable_fruit(:))    * (litterfall_fruit / lctlib_cn_leaf / lctlib_np_leaf)     &
      &                                 + (1.0_wp - fp_soluable_seed_bed(:)) * (litterfall_seedbed / lctlib_cn_leaf / lctlib_np_leaf)) &
      &                              / soil_depth_sl(:,1)
    ! sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) = &
      ! &                              sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) &
      ! &                              + ((1.0_wp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC13, :)      &
      ! &                                 + (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC13, :)     &
      ! &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixC13, :)) &
      ! &                              / soil_depth_sl(:,1)
    ! sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) = &
      ! &                              sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) &
      ! &                              + ((1.0_wp - fc_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixC14, :)      &
      ! &                                 + (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixC14, :)     &
      ! &                                 + (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixC14, :)) &
      ! &                              / soil_depth_sl(:,1)
    ! sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) = &
      ! &                              sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) &
      ! &                              + ((1.0_wp - fn_soluable_leaf(:))       * veg_litterfall_mt(ix_leaf, ixN15, :)      &
      ! &                                 + (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall_mt(ix_fruit, ixN15, :)     &
      ! &                                 + (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall_mt(ix_seed_bed, ixN15, :)) &
      ! &                              / soil_depth_sl(:,1)

    !>4.0 Add surface litter fall from sapwood and heartwood pools to surface woody litter pool for trees
    !>    for grasses, stems and halms (=sapwood) are not actually wood, and therefore paritioned to
    !>    metabolic and structural litter as if their were fine roots
    SELECT CASE(lctlib_growthform)
    ! tree PFT
    CASE (ITREE)
      ! woody litter
      ! @TODO may be re-written to a one-liner using bgcm matrices
      sb_formation_mt(ix_woody_litter, ixC, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixC, :, 1) &
        &           + litterfall_wood / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixN, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixN, :, 1) &
        &           + (litterfall_wood / wood_litter_cn) / soil_depth_sl(:,1)
      sb_formation_mt(ix_woody_litter, ixP, :, 1) = &
        &           sb_formation_mt(ix_woody_litter, ixP, :, 1) &
        &           + (litterfall_wood / wood_litter_cn / wood_litter_np) / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_woody_litter, ixC13, :, 1) = &
        ! &           sb_formation_mt(ix_woody_litter, ixC13, :, 1) &
        ! &           + (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_woody_litter, ixC14, :, 1) = &
        ! &           sb_formation_mt(ix_woody_litter, ixC14, :, 1) &
        ! &           + (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_woody_litter, ixN15, :, 1) = &
        ! &           sb_formation_mt(ix_woody_litter, ixN15, :, 1) &
        ! &           + (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) / soil_depth_sl(:,1)
    ! grass PFT
    CASE (IGRASS,ICROP)
      ! soluable litter
      sb_formation_mt(ix_soluable_litter, ixC, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixC, :, 1) &
        &             + fc_soluable_fine_root(:) &
        &             * litterfall_wood &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixN, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixN, :, 1) &
        &             + fn_soluable_fine_root(:) &
        &             * (litterfall_wood / wood_litter_cn) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_soluable_litter, ixP, :, 1) = &
        &             sb_formation_mt(ix_soluable_litter, ixP, :, 1) &
        &             + fp_soluable_fine_root(:) &
        &             * (litterfall_wood / wood_litter_cn / wood_litter_np) &
        &             / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_soluable_litter, ixC13, :, 1) = &
        ! &             sb_formation_mt(ix_soluable_litter, ixC13, :, 1) &
        ! &             + fc_soluable_fine_root(:) &
        ! &             * (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) &
        ! &             / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_soluable_litter, ixC14, :, 1) = &
        ! &             sb_formation_mt(ix_soluable_litter, ixC14, :, 1) &
        ! &             + fc_soluable_fine_root(:) &
        ! &             * (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) &
        ! &             / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_soluable_litter, ixN15, :, 1) = &
        ! &             sb_formation_mt(ix_soluable_litter, ixN15, :, 1) &
        ! &             + fn_soluable_fine_root(:) &
        ! &             * (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) &
        ! &             / soil_depth_sl(:,1)
      ! polymeric litter
      sb_formation_mt(ix_polymeric_litter, ixC, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixC, :, 1) &
        &             + (1.0_wp - fc_soluable_fine_root(:)) &
        &             * litterfall_wood &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixN, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixN, :, 1) &
        &             + (1.0_wp - fn_soluable_fine_root(:)) &
        &             * (litterfall_wood / wood_litter_cn) &
        &             / soil_depth_sl(:,1)
      sb_formation_mt(ix_polymeric_litter, ixP, :, 1) = &
        &             sb_formation_mt(ix_polymeric_litter, ixP, :, 1) &
        &             + (1.0_wp - fp_soluable_fine_root(:)) &
        &             * (litterfall_wood / wood_litter_cn / wood_litter_np) &
        &             / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) = &
        ! &             sb_formation_mt(ix_polymeric_litter, ixC13, :, 1) &
        ! &             + (1.0_wp - fc_soluable_fine_root(:)) &
        ! &             * (veg_litterfall_mt(ix_sap_wood, ixC13, :) + veg_litterfall_mt(ix_heart_wood, ixC13, :)) &
        ! &             / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) = &
        ! &             sb_formation_mt(ix_polymeric_litter, ixC14, :, 1) &
        ! &             + (1.0_wp - fc_soluable_fine_root(:)) &
        ! &             * (veg_litterfall_mt(ix_sap_wood, ixC14, :) + veg_litterfall_mt(ix_heart_wood, ixC14, :)) &
        ! &             / soil_depth_sl(:,1)
      ! sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) = &
        ! &             sb_formation_mt(ix_polymeric_litter, ixN15, :, 1) &
        ! &             + (1.0_wp - fn_soluable_fine_root(:)) &
        ! &             * (veg_litterfall_mt(ix_sap_wood, ixN15, :) + veg_litterfall_mt(ix_heart_wood, ixN15, :)) &
        ! &             / soil_depth_sl(:,1)
    END SELECT

    !>5.0 Distribute fine and coarse root litter to soil depth according to the root profile
    !> litterfall_fineroot, litterfall_coarseroot
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        ! soluable litter (fine root & coarse root)
        sb_formation_mt(ix_soluable_litter, ixC, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixC, ic, is) &
          &             + (fc_soluable_fine_root(ic)   * litterfall_fineroot &
          &             +  fc_soluable_coarse_root(ic) * litterfall_coarseroot) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixN, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixN, ic, is) &
          &             + (fn_soluable_fine_root(ic)   * (litterfall_fineroot/fineroot_litter_cn) &
          &             +  fn_soluable_coarse_root(ic) * (litterfall_coarseroot/coarseroot_litter_cn)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        sb_formation_mt(ix_soluable_litter, ixP, ic, is) = &
          &             sb_formation_mt(ix_soluable_litter, ixP, ic, is) &
          &             + (fp_soluable_fine_root(ic)   * (litterfall_fineroot/fineroot_litter_cn/fineroot_litter_np) &
          &             +  fp_soluable_coarse_root(ic) * (litterfall_coarseroot/coarseroot_litter_cn/coarseroot_litter_np)) &
          &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        ! sb_formation_mt(ix_soluable_litter, ixC13, ic, is) = &
          ! &             sb_formation_mt(ix_soluable_litter, ixC13, ic, is) &
          ! &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC13, ic) &
          ! &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
          ! &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        ! sb_formation_mt(ix_soluable_litter, ixC14, ic, is) = &
          ! &             sb_formation_mt(ix_soluable_litter, ixC14, ic, is) &
          ! &             + (fc_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixC14, ic) &
          ! &             +  fc_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
          ! &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        ! sb_formation_mt(ix_soluable_litter, ixN15, ic, is) = &
          ! &             sb_formation_mt(ix_soluable_litter, ixN15, ic, is) &
          ! &             + (fn_soluable_fine_root(ic)   * veg_litterfall_mt(ix_fine_root, ixN15, ic) &
          ! &             +  fn_soluable_coarse_root(ic) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
          ! &             * root_fraction_sl(ic,is)  / soil_depth_sl(ic,is)
        ! simple soil model / jsm
        SELECT CASE(TRIM(sb_model_scheme))
        CASE ("simple_1d")
          ! polymeric litter (fine root & coarse root)
          sb_formation_mt(ix_polymeric_litter, ixC, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic))   * litterfall_fineroot &
            &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * litterfall_coarseroot) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic))   * (litterfall_fineroot/fineroot_litter_cn) &
            &           +  (1.0_wp - fn_soluable_coarse_root(ic)) * (litterfall_coarseroot/coarseroot_litter_cn)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixP, ic, is)    = &
            &           sb_formation_mt(ix_polymeric_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_fine_root(ic))   * (litterfall_fineroot/fineroot_litter_cn/fineroot_litter_np) &
            &           +  (1.0_wp - fp_soluable_coarse_root(ic)) * (litterfall_coarseroot/coarseroot_litter_cn/coarseroot_litter_np)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_polymeric_litter, ixC13, ic, is)    = &
            ! &           sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) &
            ! &           + ((1.0_wp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC13, ic) &
            ! &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_polymeric_litter, ixC14, ic, is)    = &
            ! &           sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) &
            ! &           + ((1.0_wp - fc_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixC14, ic) &
            ! &           +  (1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_polymeric_litter, ixN15, ic, is)    = &
            ! &           sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) &
            ! &           + ((1.0_wp - fn_soluable_fine_root(ic))   * veg_litterfall_mt(ix_fine_root, ixN15, ic) &
            ! &           +  (1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
        CASE ("jsm")
          ! polymeric litter (fine root)
          sb_formation_mt(ix_polymeric_litter, ixC, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_fine_root(ic)) * litterfall_fineroot) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixN, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_fine_root(ic)) * (litterfall_fineroot/fineroot_litter_cn)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_polymeric_litter, ixP, ic, is) = &
            &           sb_formation_mt(ix_polymeric_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_fine_root(ic)) * (litterfall_fineroot/fineroot_litter_cn/fineroot_litter_np)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) = &
            ! &           sb_formation_mt(ix_polymeric_litter, ixC13, ic, is) &
            ! &           + ((1.0_wp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC13, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) = &
            ! &           sb_formation_mt(ix_polymeric_litter, ixC14, ic, is) &
            ! &           + ((1.0_wp - fc_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixC14, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) = &
            ! &           sb_formation_mt(ix_polymeric_litter, ixN15, ic, is) &
            ! &           + ((1.0_wp - fn_soluable_fine_root(ic)) * veg_litterfall_mt(ix_fine_root, ixN15, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! polymeric litter (coarse root)
          sb_formation_mt(ix_woody_litter, ixC, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixC, ic, is) &
            &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * litterfall_coarseroot) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixN, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixN, ic, is) &
            &           + ((1.0_wp - fn_soluable_coarse_root(ic)) * (litterfall_coarseroot/coarseroot_litter_cn)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          sb_formation_mt(ix_woody_litter, ixP, ic, is) = &
            &           sb_formation_mt(ix_woody_litter, ixP, ic, is) &
            &           + ((1.0_wp - fp_soluable_coarse_root(ic)) * (litterfall_coarseroot/coarseroot_litter_cn/coarseroot_litter_np)) &
            &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_woody_litter, ixC13, ic, is) = &
            ! &           sb_formation_mt(ix_woody_litter, ixC13, ic, is) &
            ! &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC13, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_woody_litter, ixC14, ic, is) = &
            ! &           sb_formation_mt(ix_woody_litter, ixC14, ic, is) &
            ! &           + ((1.0_wp - fc_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixC14, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
          ! sb_formation_mt(ix_woody_litter, ixN15, ic, is) = &
            ! &           sb_formation_mt(ix_woody_litter, ixN15, ic, is) &
            ! &           + ((1.0_wp - fn_soluable_coarse_root(ic)) * veg_litterfall_mt(ix_coarse_root, ixN15, ic)) &
            ! &           * root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
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

  ! --- convert_lct_to_sla() ---------------------------------------------------------------------------------------------------
  !
  ! find the sla value for JSBACH lctNumber
  !  
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE quincy_allometry_turnover(  &
         &                    pheno_type, &
         &                    C_2_litter_leaf, C_2_litter_wood, &
         &                    c_green, c_woods, c_reserve, &
         &                    Cflx_NPP_2_green_pot, Cflx_NPP_2_wood_pot, Cflx_NPP_2_reserve_pot, &
         &                    c_leaf, c_fineroot, c_fruit, c_seedbed, c_sapwood, c_heartwood, c_coarseroot, &
         &                    lit_leaf, lit_fineroot, &
         &                    lit_fruit, lit_seedbed, &
         &                    lit_wood, lit_coarseroot)

    USE mo_veg_constants,           ONLY: sm2lm_grass, background_mort_rate_tree, background_mort_rate_grass, 
    IMPLICIT NONE
    INTEGER,  INTENT(in)    :: pheno_type             ! crop =5
    !
    REAL(wp), INTENT(inout) :: c_green            ! Green carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_woods            ! Wood carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_reserve          ! Reserve carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_leaf             ! Leaf carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_fineroot         ! Fineroot carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_fruit            ! Fruit carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_seedbed          ! Seedbed carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_sapwood          ! SapWood carbon pool (sap and heart wood combined): on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_heartwood        ! HeartWood carbon pool (sap and heart wood combined): on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_coarseroot       ! Coarseroot carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(in) :: C_2_litter_leaf       ! litter flux to leaf and fineroot [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in) :: C_2_litter_wood       ! litter flux to wood and coarse root [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: Cflx_NPP_2_green_pot  ! carbon flux to leaf and fineroot [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: Cflx_NPP_2_wood_pot   ! carbon flux to reserve pool [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: Cflx_NPP_2_reserve_pot ! carbon flux to wood and coarse root [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: lit_leaf       ! carbon flux to leaf and fineroot [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: lit_fineroot   ! carbon flux to reserve pool [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: lit_fruit      ! carbon flux to wood and coarse root [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: lit_seedbed    ! carbon flux to leaf and fineroot [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: lit_wood       ! carbon flux to reserve pool [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: lit_coarseroot ! carbon flux to wood and coarse root [mol(C)/m^2(canopy)/yr]


    ! locals
    REAL(wp) :: lit_leaf_pot, lit_fineroot_pot, lit_fruit_pot, lit_seedbed_pot, &
                lit_coarseroot_pot, lit_sapwood_pot, lit_heartwood_pot
    REAL(wp) :: jsb_allocfrac_leaf, jsb_allocfrac_fineroot, jsb_allocfrac_fruit, &
                jsb_allocfrac_coarseroot, jsb_allocfrac_wood
    REAL(wp) :: tau_leaf, tau_fineroot, tau_seedbed, tau_coarseroot, tau_branch, tau_sapwood, &
                frac_sapwood_branch, mortality
    REAL(wp) :: hlp1

    !> init local variables
    tau_leaf       = 0._wp
    tau_fineroot   = 0._wp
    tau_seedbed    = 0._wp
    tau_coarseroot = 0._wp
    tau_branch     = 0._wp
    tau_sapwood    = 0._wp
    mortality      = 0._wp
    frac_sapwood_branch = 0._wp
    
    !> the litter is 
    SELECT CASE(PhenologyType)
    CASE(1)
    ! evergreen trees
      jsb_allocfrac_leaf     =
      jsb_allocfrac_fineroot =
      jsb_allocfrac_fruit    =
      jsb_litfrac_coarseroot = 
      mortality              = background_mort_rate_tree
    CASE(2)
    ! summergreen trees
      jsb_litfrac_leaf       =
      jsb_litfrac_fruit      =
      jsb_litfrac_fineroot   =
      jsb_litfrac_coarseroot = 
      mortality              = background_mort_rate_tree
    CASE(3)
    ! raingreen trees
      jsb_litfrac_leaf       =
      jsb_litfrac_fruit      =
      jsb_litfrac_fineroot   =
      jsb_litfrac_coarseroot = 
      mortality              = background_mort_rate_tree
    CASE(4)
    ! grasses
      jsb_litfrac_leaf       =
      jsb_litfrac_fruit      =
      jsb_litfrac_fineroot   =
      jsb_litfrac_coarseroot = 
      mortality              = background_mort_rate_grass
    CASE(5)
    ! crops
      jsb_litfrac_leaf       =
      jsb_litfrac_fruit      =
      jsb_litfrac_fineroot   =
      jsb_litfrac_coarseroot = 
      mortality              = 0._wp ! crops are harvested every year
    END SELECT
  


    !! ============== UPDATE ALL LITTER FLUXES =================================================================
    ! potential litter flux is calculated by the pool size, turnover rate, and mortality (QUINCY defined)
    lit_leaf_pot     = c_leaf     * (1._wp/tau_leaf + mortality)
    lit_fineroot_pot = c_fineroot * (1._wp/tau_fineroot + mortality)
    lit_fruit_pot    = c_fruit    * mortality
    lit_seedbed_pot  = c_seedbed  * (1._wp/tau_seedbed + mortality)
    
    lit_coarseroot_pot = c_coarseroot   * (1._wp/tau_coarseroot + mortality)
    lit_sapwood_pot    = c_sapwood * (1._wp/tau_branch*frac_sapwood_branch + mortality)
    lit_heartwood_pot  = c_heartwood * mortality
    
    ! actual litter flux is constrained by C_2_litter_leaf and C_2_litter_wood
    ! C_2_litter_leaf
    IF (C_2_litter_leaf > eps8) THEN
      IF ( (lit_leaf_pot + lit_fineroot_pot + lit_fruit_pot + lit_seedbed_pot) > C_2_litter_leaf ) THEN
      
      END If
    ELSE
      lit_leaf     = 0._wp
      lit_fineroot = 0._wp
      lit_fruit    = 0._wp
      lit_seedbed  = 0._wp
    END If
    ! C_2_litter_wood
    IF (C_2_litter_wood > eps8) THEN
      IF ( (lit_coarseroot_pot + lit_sapwood_pot + lit_heartwood_pot) > C_2_litter_wood ) THEN
        
      END IF
    ELSE
      lit_sapwood    = 0._wp
      lit_heartwood  = 0._wp
      lit_coarseroot = 0._wp
    END If   
    
    lit_fruit      = C_2_litter_leaf * jsb_litfrac_fruit
    lit_fineroot   = C_2_litter_leaf * jsb_litfrac_fineroot
    lit_leaf       = C_2_litter_leaf * jsb_litfrac_leaf
    lit_seedbed    = c_seedbed * tau_seedbed
    
    lit_coarseroot = C_2_litter_wood * jsb_litfrac_coarseroot
    lit_wood       = C_2_litter_wood - lit_coarseroot  

    !! ============== UPDATE ALL C-POOLS =======================================================================
    ! npp to green pools QUINCY
    c_leaf         = c_leaf         + Cflx_NPP_2_green_pot * jsb_allocfrac_leaf     - lit_leaf
    c_fineroot     = c_fineroot     + Cflx_NPP_2_green_pot * jsb_allocfrac_fineroot - lit_fineroot
    c_fruit        = c_fruit        + Cflx_NPP_2_green_pot * jsb_allocfrac_fruit    - lit_fruit - c_fruit * tau_fruit
    c_seedbed      = c_seedbed      + c_fruit * tau_fruit                           - lit_seedbed
    
    ! npp to wood pools QUINCY
    c_coarseroot   = c_coarseroot   + Cflx_NPP_2_wood_pot  * jsb_allocfrac_coarseroot - lit_coarseroot
    c_sapwood      = c_sapwood      + Cflx_NPP_2_wood_pot  * jsb_allocfrac_wood       - lit_sapwood - &
            &         (1._wp - frac_sapwood_branch) * tau_sapwood
    c_heartwood    = c_heartwood    + (1._wp - frac_sapwood_branch) * tau_sapwood     - lit_heartwood
    
    ! JSBACH pools
    c_green       = c_green       + Cflx_NPP_2_green_pot
    c_woods       = c_woods       + Cflx_NPP_2_wood_pot
    c_reserve     = c_reserve     + Cflx_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N
    

  END SUBROUTINE quincy_allometry_turnover

  ! --- convert_lct_to_sla() ---------------------------------------------------------------------------------------------------
  !
  ! find the sla value for JSBACH lctNumber
  !  
  SUBROUTINE convert_jsb_lct_to_quincy(jsb_LctNumber, lctlib_sla, lctlib_cn_leaf, lctlib_np_leaf, &
                 &                     lctlib_tau_leaf, lctlib_tau_fine_root, lctlib_tau_coarse_root, &
                 &                     lctlib_tau_branch, lctlib_tau_wood, lctlib_tau_fruit, lctlib_tau_seedbed, &
                 &                     lctlib_k_crtos)
  ! #  QUINCY LCT Number                            
  ! #  1: BEM:  moist broadleaved evergreen
  ! #  2: BED:  dry broadleaved evergreen
  ! #  3: BDR:  rain green broadleaved deciduous
  ! #  4: BDS:  summer green broadleaved deciduous
  ! #  5: NE :  needle-leaved evergreen
  ! #  6: NS :  summer green needle-leaved
  ! #  7: TeH:  C3 grass
  ! #  8: TrH:  C4 grass
  ! #  9: TeP:  C3 pasture
  ! # 10: TrP:  C4 pasture
  ! # 11: TeC:  C3 crop
  ! # 12: TrC:  C4 crop
  ! # 13: BSO:  bare soil
  ! # 14: UAR:  urban area  
  !------------------------------------------------------------------------------------
  ! QUINCY sla
  ! # From Kattge et al. 2011, Table 5, converted internally from (mm/mgDW) to (m2/mol C)        
  ! sla 9.46 8.99 15.4 15.4 5.0 10.09 20.12 19.23 20.12 19.23 20.12 19.23 19.23 19.23
  !------------------------------------------------------------------------------------
  ! QUINCY cn_leaf: default leaf nitrogen concentration (mgN/gDW)        
  ! # From Kattge et al. 2011, Table 5, converted internally from nitrogen concentration (mg/gDW) to CN mol C / mol N        
  ! cn_leaf 16.89 13.73 21.32 21.32 12.09 19.37 17.84 14.14 17.84 14.14 17.84 14.14 14.14 14.14
  !------------------------------------------------------------------------------------
  ! QUINCY k_crtos: coarse root to sapwood mass ratio         
  !k_crtos 0.1 0.1 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2
  !------------------------------------------------------------------------------------
  ! # JSBACH LCT Number
  ! #  1: Glacier                     
  ! #  2: Tropical evergreen trees                     
  ! #  3: Tropical deciduous trees                     
  ! #  4: Extra-tropical evergreen trees                     
  ! #  5: Extra-tropical deciduous trees                     
  ! #  6: Temperate broadleaf evergreen trees                     
  ! #  7: Temperate broadleaf deciduous trees                     
  ! #  8: Coniferous evergreen trees                     
  ! #  9: Coniferous deciduous trees                     
  ! # 10: Raingreen shrubs                     
  ! # 11: Deciduous shrubs                     
  ! # 12: C3 grass                     
  ! # 13: C4 grass                     
  ! # 14: Pasture                     
  ! # 15: C3 pasture                     
  ! # 16: C4 pasture                     
  ! # 17: Tundra                     
  ! # 18: Swamp                     
  ! # 19: Crops                     
  ! # 20: C3 crop                     
  ! # 21: C4 crop   
  !--------------------------------------------------
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER, INTENT(in)                       :: jsb_LctNumber            !< lct number from JSBACH lctlib file
    INTEGER, INTENT(out)                      :: quincy_LctNumber         !< corresponding QUINCY lct number
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_sla               !< corresponding QUINCY lct sla
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_cn_leaf           !< corresponding QUINCY lct cn_leaf
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_np_leaf           !< corresponding QUINCY lct np_leaf
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_k_crtos           !< corresponding QUINCY lct k_crtos
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_leaf          !< corresponding QUINCY lct tau_leaf
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_fine_root     !< corresponding QUINCY lct tau_fine_root
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_coarse_root   !< corresponding QUINCY lct tau_coarse_root
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_branch        !< corresponding QUINCY lct tau_branch
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_wood          !< corresponding QUINCY lct tau_sap_wood
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_fruit         !< corresponding QUINCY lct tau_fruit
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_seedbed       !< corresponding QUINCY lct tau_seed_litter
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_k_crtos           !< corresponding QUINCY lct k_crtos
    ! ---------------------------
    ! 0.2 Local

    !> init out var
    !>
    
    SELECT CASE(jsb_LctNumber)
    CASE(2, 6)
      !> moist broadleaved evergreen
      quincy_LctNumber = 1
      IF (PRESENT(lctlib_sla)) lctlib_sla     = 9.46_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(4)
      !> dry broadleaved evergreen
      quincy_LctNumber = 2
      IF (PRESENT(lctlib_sla)) lctlib_sla     = 8.99_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(3, 5, 10)
      !> rain green broadleaved deciduous
      quincy_LctNumber = 3
      IF (PRESENT(lctlib_sla)) lctlib_sla     = 15.4_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(7, 11)
      !> summer green broadleaved deciduous
      quincy_LctNumber = 4
      IF (PRESENT(lctlib_sla)) lctlib_sla     = 15.4_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(8)
      !> needle-leaved evergreen
      quincy_LctNumber = 5
      IF (PRESENT(lctlib_sla)) lctlib_sla     = 5.0_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(9)
      !> summer green needle-leaved
      quincy_LctNumber = 6
      IF (PRESENT(lctlib_sla)) lctlib_sla     = 10.09_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(12)
      !> C3 grass
      quincy_LctNumber = 7
      IF (PRESENT(lctlib_sla)) lctlib_sla = 20.12_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(14, 15)
      !> C3 pasture
      quincy_LctNumber = 9
      IF (PRESENT(lctlib_sla)) lctlib_sla = 20.12_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(19, 20)
      !> C3 crop
      quincy_LctNumber = 11
      IF (PRESENT(lctlib_sla)) lctlib_sla = 20.12_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(13)
      !> C4 grass
      quincy_LctNumber = 8
      IF (PRESENT(lctlib_sla)) lctlib_sla = 20.12_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(16)
      !> C4 pasture
      quincy_LctNumber = 10
      IF (PRESENT(lctlib_sla)) lctlib_sla = 20.12_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(21)
      !> C4 crop
      quincy_LctNumber = 12
      IF (PRESENT(lctlib_sla)) lctlib_sla = 20.12_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    CASE(1, 17, 18)
      !> glacier, Tundra, Swamp
      quincy_LctNumber = 13
      IF (PRESENT(lctlib_sla)) lctlib_sla = 19.23_wp
      IF (PRESENT(lctlib_cn_leaf)) lctlib_cn_leaf = 
      IF (PRESENT(lctlib_k_crtos)) lctlib_k_crtos = 
    END SELECT          
  
  END SUBROUTINE convert_jsb_lct_to_quincy
  
  
  SUBROUTINE calc_litter_fraction ( &
    &  alloc_fract1, tau1, & 
    &  alloc_fract2, tau2, & 
    &  mature_year, &
    &  lit_fract1,  &
    &  lit_fract2,  &
    &  alloc_fract3, tau3, &     
    &  lit_fract3)
    
    USE mo_kind,                        ONLY: wp
    USE mo_jsb_math_constants,          ONLY: eps8
    USE mo_veg_constants,               ONLY: k_f_demand, omega_nutrient_demand
    
    IMPLICIT NONE
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp),               INTENT(in)  :: alloc_fract1        !< allocation fraction of organism one
    REAL(wp),               INTENT(in)  :: alloc_fract2        !< allocation fraction of organism two
    REAL(wp),               INTENT(in)  :: tau1                !< turnover rate of organism one [y -1]
    REAL(wp),               INTENT(in)  :: tau2                !< turnover rate of organism one [y -2]
    REAL(wp),               INTENT(in)  :: mature_year         !< number of years for the PhenologyType to mature
    REAL(wp),               INTENT(inout)  :: lit_fract1       !< fraction of litter from organism one
    REAL(wp),               INTENT(inout)  :: lit_fract2       !< fraction of litter from organism two
    REAL(wp), OPTIONAL      INTENT(in)  :: alloc_fract3        !< allocation fraction of organism three
    REAL(wp), OPTIONAL      INTENT(in)  :: tau3                !< turnover rate of organism three [y -1]
    REAL(wp), OPTIONAL      INTENT(inout)  :: lit_fract3       !< fraction of litter from organism three
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: org_pool1,org_pool2,org_pool3
    INTEGER,                    :: yr
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_fraction'
    
    !> init local variables
    org_pool1 = 0._wp
    org_pool2 = 0._wp
    org_pool3 = 0._wp
    
    !> 1.0 
    DO yr = 1, mature_year
      org_pool1 = (org_pool1 + alloc_fract1) * (1._wp - 1._wp/tau1)
      org_pool2 = (org_pool2 + alloc_fract2) * (1._wp - 1._wp/tau2)
      IF PRESENT(alloc_fract3) THEN
        org_pool3 = (org_pool3 + alloc_fract3) * (1._wp - 1._wp/tau3)
      END If
    END DO

  
  END SUBROUTINE calc_litter_fraction
  
END MODULE mo_sb_carbon_process
