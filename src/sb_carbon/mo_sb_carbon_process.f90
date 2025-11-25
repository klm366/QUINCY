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

  PUBLIC :: calc_jsb_C_litter, prep_jsb_sb_processes, get_per_tile, convert_lct_to_sla, calc_sb_treeheight

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_carbon_process'

CONTAINS

  ! Calculations to allocate NPP to the plant C pools
  ! Step 1: allocate NPP to SB litter pools
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_jsb_C_litter( &
    & jsb_allocation,           &,! in
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
    & jsb_lctnumber,              & ! in
    !
    & c_green,                    & ! inout
    & c_woods,                    & ! inout
    & c_reserve,                  & ! inout
    & c_crop_harvest,             & ! inout
    & c_leaf,                     & ! inout
    & c_fineroot,                 & ! inout
    & c_fruit,                    & ! inout
    & c_seedbed,                  & ! inout
    & c_sapwood,                  & ! inout
    & c_heartwood,                & ! inout
    & c_coarseroot,               & ! inout
    & tree_height,                & ! inout
    & f_lai,                      & ! inout
    & dheight,                    & ! inout
    & df_lai,                     & ! inout
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
    USE mo_veg_constants,          ONLY: background_mort_rate_tree, background_mort_rate_grass
    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    ! Use of process memory
    dsl4jsb_Use_memory(SB_)    

    ! Arguments
    LOGICAL,  INTENT(in)    :: jsb_allocation         ! use jsb allocation with fixed allocation from quincy
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
    INTEGER,  INTENT(in)    :: jsb_LctNumber          ! lct number in JSBACH
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
    REAL(wp), INTENT(inout) :: tree_height       ! tree height [m]
    REAL(wp), INTENT(inout) :: f_lai             ! lai index to estimate tree height
    REAL(wp), INTENT(inout) :: dheight       ! tree height [m]
    REAL(wp), INTENT(inout) :: df_lai             ! lai index to estimate tree height
    
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
    REAL(wp), INTENT(inout)   :: litterfall_sapwood            ! amount of leafLitter
    REAL(wp), INTENT(inout)   :: litterfall_heartwood            ! amount of leafLitter
    REAL(wp), INTENT(inout)   :: litterfall_coarseroot   ! amount of woodyLitter

    !
    ! added for thonicke fire algorithm (spitfire) needed to correctly split the carbon pools into the different fuel classes
    !real(wp),optional,intent(inout) :: fract_litter_wood_new   ! new fraction in above ground wood litter pool


    ! locals
    LOGICAL  :: is_crop                 ! logical crop mask
    REAL(wp) :: c_green_max,c_reserve_optimal
    REAL(wp) :: Cflx_2_crop_harvest, Cflx_crop_harvest_2_atmos
    REAL(wp) :: Cflx_faeces_2_LG, Cflx_faeces_2_atm
    REAL(wp) :: excess_carbon
    ! R: only for spitfire
    !REAL(wp) :: litter_wood_new, litter_wood_total

    ! local variables
    REAL(wp) :: NPP_2_green_via_mobileN, NPP_2_wood_via_mobileN
    REAL(wp) :: redFact_Nlimitation, red_scalar
    
    REAL(wp) :: lit_leaf_pot, lit_fineroot_pot, lit_fruit_pot, lit_seedbed_pot, &
                lit_coarseroot_pot, lit_sapwood_pot, lit_heartwood_pot
    REAL(wp) :: q_tau_leaf, q_tau_fineroot, q_tau_seedbed, q_tau_coarseroot, q_tau_branch, q_tau_sapwood, q_tau_fruit, &
                q_mortality, frac_sapwood_branch


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

    NPP_act                  = 0.0_wp

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

    q_tau_leaf       = 0._wp
    q_tau_fineroot   = 0._wp
    q_tau_seedbed    = 0._wp
    q_tau_coarseroot = 0._wp
    q_tau_branch     = 0._wp
    q_tau_sapwood    = 0._wp
    q_tau_fruit      = 0._wp
    q_mortality      = 0._wp
    frac_sapwood_branch = 0._wp
    
    !> retrieve QUINCY lctlib values for turnover time
    CALL convert_jsb_lct_to_quincy(jsb_LctNumber, quincy_LctNumber, &
      &   lctlib_tau_leaf= q_tau_leaf, &
      &   lctlib_tau_fine_root = q_tau_fineroot , &
      &   lctlib_tau_seedbed = q_tau_seedbed, &
      &   lctlib_tau_coarse_root = q_tau_coarseroot, &
      &   lctlib_tau_branch = q_tau_branch, &
      &   lctlib_tau_wood = q_tau_sapwood, &
      &   lctlib_tau_fruit = q_tau_fruit &
      &   lctlib_frac_sapwood_branch = frac_sapwood_branch &
      &   )
    SELECT CASE(pheno_type)
    CASE(1, 2, 3)
    ! trees
      q_mortality = background_mort_rate_tree
    CASE(4)
    ! grasses
      q_mortality = background_mort_rate_grass
    CASE(5)
    ! crops
      q_mortality = 0._wp ! crops are harvested every year
    END SELECT
    ! on all veg tiles:

    c_green_max       = greenC2leafC * LAI / specific_leaf_area_C
    c_reserve_optimal = reserveC2leafC * max_LAI / specific_leaf_area_C

    !! ============== NON-N-LIMITED CARBON ALLOCATION =========================================================
    !!
    !! Perform those C-allocation steps that are not restricted by N-availability
    !! 1. Leaf shedding: follow the jsb assumption for green pool shedding
    !! 2. Herbivory loss: jsb or quincy
    !! 3. Wood shedding: jsb or quincy
    !! 4. Depletion of the reserve pool: jsb
    !!
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
    !! 2. Herbivory loss from c_green & Npool_green w.r.t PFTs in lctlib

    IF (jsb_allocation) THEN
      ! JSBACH
      Green_2_herbivory  = c_green * fract_green_2_herbivory             !! Grazing flux indirect comput.from NPP losses
    ELSE
      ! QUINCY
      Green_2_herbivory  = (c_leaf + c_fruit) * fract_green_2_herbivory  !! Herbivory loss only from leaf and fruit
      c_leaf = c_leaf * (1._wp - fract_green_2_herbivory)
      c_fruit = c_fruit * (1._wp - fract_green_2_herbivory)
    END IF
    c_green        = c_green - Green_2_herbivory                       !! Loss due to grazing

    Cflx_faeces_2_LG  = Green_2_herbivory * fract_C_faeces2_LG             !! Cflux_faeces_2_LG to be stored in litter green pool
                                                                           !!      (not considered in cflux_c_greenwood_2_litter)
                                                                           !!      (no life times (tau) enter Cflx_faeces_2_LG!)
    Cflx_faeces_2_atm = Green_2_herbivory * (1.0_wp - fract_C_faeces2_LG)  !! Cflux_faeces_2_atm to be lost to atmosphere
                                                                           !!      (dk: not added to CO2?)
    !! ---------------------------------------------------------------------------------
    !! 3. Wood shedding: Transfer of C from wood pool to wood litter pools
    !!
    
    ! QUINCY veg pool potential litter flux
    lit_coarseroot_pot = c_coarseroot * (1._wp/q_tau_coarseroot + q_mortality)
    lit_sapwood_pot    = c_sapwood    * (1._wp/q_tau_branch*frac_sapwood_branch + q_mortality)
    lit_heartwood_pot  = c_heartwood  * q_mortality    
    IF (jsb_allocation) THEN
      ! JSBACH
                                                              !! Assuming that forests continuously die at the inverse lifetime of trees
      C_2_litter_woodPool = c_woods / tau_c_woods             ! The bigger tau gets, the weaker the C pools will be depreciated:
                                                              ! tau_c_woods: [days] -> [timestep]
                                                              !! .. the shedded wood (MAX FUNCTION FOR NONWOODY PFTS)
      !> C_2_litter_wood
      IF (C_2_litter_woodPool > eps8) THEN
        IF ( (lit_coarseroot_pot + lit_sapwood_pot + lit_heartwood_pot) > C_2_litter_woodPool ) THEN
          red_scalar = C_2_litter_woodPool / (lit_coarseroot_pot + lit_sapwood_pot + lit_heartwood_pot)
          litterfall_sapwood    = lit_sapwood_pot    * red_scalar
          litterfall_heartwood  = lit_heartwood_pot  * red_scalar
          litterfall_coarseroot = lit_coarseroot_pot * red_scalar
        END IF
      ELSE
        litterfall_sapwood    = 0._wp
        litterfall_heartwood  = 0._wp
        litterfall_coarseroot = 0._wp
      END If
      ! R: only for spitfire
      ! if (present(fract_litter_wood_new))     litter_wood_new = fract_wood_aboveGround * C_2_litter_woodPool
    ELSE
      ! QUINCY
      litterfall_sapwood    = 0._wp
      litterfall_heartwood  = 0._wp
      litterfall_coarseroot = 0._wp
      C_2_litter_woodPool = litterfall_sapwood + litterfall_heartwood + litterfall_coarseroot
    END IF
    c_woods      = c_woods -  C_2_litter_woodPool                 !! .. and then all is subtracted from the wood pool.

    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter + C_2_litter_woodPool

    !! ----------------------------------------------------------------------------
    !! 4. Depletion of reserve pool: Transfer of C from reserve pool to green litter pool
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

    leafLitter    = cflux_c_greenwood_2_litter - C_2_litter_woodPool + Cflx_faeces_2_LG - Cflx_2_crop_harvest
    
    
    ! QUINCY distrubution for green pool litter
    !> potential litter flux is calculated by the pool size, turnover rate, and mortality (QUINCY defined)
    lit_leaf_pot     = c_leaf     * (1._wp/q_tau_leaf + q_mortality)
    lit_fineroot_pot = c_fineroot * (1._wp/q_tau_fineroot + q_mortality)
    lit_fruit_pot    = c_fruit    * q_mortality
    lit_seedbed_pot  = c_seedbed  * (1._wp/q_tau_seedbed + q_mortality)
    
    !> actual litter flux is constrained by leafLitter
    IF (leafLitter > eps8) THEN
      IF ( (lit_leaf_pot + lit_fineroot_pot + lit_fruit_pot + lit_seedbed_pot) > leafLitter ) THEN
        red_scalar = leafLitter / (lit_leaf_pot + lit_fineroot_pot + lit_fruit_pot + lit_seedbed_pot)
        litterfall_leaf     = lit_leaf_pot    * red_scalar
        litterfall_fineroot = lit_fineroot    * red_scalar
        litterfall_fruit    = lit_fruit_pot   * red_scalar
        litterfall_seedbed  = lit_seedbed_pot * red_scalar
      END If
    ELSE
      litterfall_leaf     = 0._wp
      litterfall_fineroot = 0._wp
      litterfall_fruit    = 0._wp
      litterfall_seedbed  = 0._wp
    END If  
    
    !! ============== PARTITION JSB C & LITTER FLUX TO QUINCY ===================================================
    CALL sb_carbon_allometry(jsb_allocation, &
         &                    pheno_type, jsb_LctNumber, &
         &                    c_green, c_woods, c_reserve, &
         &                    NPP_pot_yDayMean * sec_per_day, NPP_act, &
         &                    fract_npp_2_woodPool, fract_npp_2_reservePool, &
         &                    fract_npp_2_exudates, Max_C_content_woods, &
         &                    LAI, max_LAI, tree_height, f_lai, dheight, df_lai， &
         ! &                    Cflx_NPP_2_green_pot, Cflx_NPP_2_wood_pot, Cflx_NPP_2_reserve_pot, &
         &                    c_leaf, c_fineroot, c_fruit, c_seedbed, c_sapwood, c_heartwood, c_coarseroot, & 
         &                    litterfall_leaf, litterfall_fineroot, &
         &                    litterfall_fruit, litterfall_seedbed, &
         &                    litterfall_sapwood, litterfall_heartwood, litterfall_coarseroot)


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
  SUBROUTINE sb_carbon_allometry(jsb_allocation, &
         &                    pheno_type, jsb_LctNumber, &
         &                    c_green, c_woods, c_reserve, &
         &                    NPP_pot, NPP_act, &
         &                    fract_npp_2_woodPool, fract_npp_2_reservePool, &
         &                    fract_npp_2_exudates, Max_C_content_woods, &
         &                    LAI, max_LAI, tree_height, f_lai, dheight, df_lai, &
         ! &                    Cflx_NPP_2_green_pot, Cflx_NPP_2_wood_pot, Cflx_NPP_2_reserve_pot, &
         &                    c_leaf, c_fineroot, c_fruit, c_seedbed, c_sapwood, c_heartwood, c_coarseroot, &
         &                    lit_leaf, lit_fineroot, &
         &                    lit_fruit, lit_seedbed, &
         &                    lit_sapwood, lit_heartwood, lit_coarseroot)

    USE mo_veg_constants,           ONLY: sm2lm_grass
    
    IMPLICIT NONE
    LOGICAL,  INTENT(in)    :: jsb_allocation       
    INTEGER,  INTENT(in)    :: pheno_type             ! crop =5
    INTEGER,  INTENT(in)    :: jsb_LctNumber          ! crop =5
    !
    REAL(wp), INTENT(in)    :: fract_npp_2_woodPool    ! Fraction of NPP to be put maximally into the green pool
    REAL(wp), INTENT(in)    :: fract_npp_2_reservePool ! Fraction of NPP to be put into the optimally into the reserve pool
    REAL(wp), INTENT(in)    :: fract_npp_2_exudates    ! Fraction of NPP to be put into the optimally into the root exudates    
    REAL(wp), INTENT(in)    :: Max_C_content_woods     ! Maximum carbon content of wood pool (from lctLib)   
    REAL(wp), INTENT(in)    :: LAI                     ! Yesterdays mean LAI
    REAL(wp), INTENT(in)    :: max_LAI                 ! maximum LAI
    REAL(wp), INTENT(in)    :: lit_leaf       ! carbon flux to leaf and fineroot [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in)    :: lit_fineroot   ! carbon flux to reserve pool [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in)    :: lit_fruit      ! carbon flux to wood and coarse root [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in)    :: lit_seedbed    ! carbon flux to leaf and fineroot [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in)    :: lit_sapwood       ! carbon flux to reserve pool [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in)    :: lit_heartwood       ! carbon flux to reserve pool [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in)    :: lit_coarseroot ! carbon flux to wood and coarse root [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in)    :: NPP_pot               ! potential NPP [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(inout) :: NPP_act               ! actual NPP [mol(C)/m^2(canopy)/yr]
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
    REAL(wp), INTENT(inout) :: tree_height        ! tree height, [m]
    REAL(wp), INTENT(inout) :: f_lai              ! lai index to estimate tree height
    REAL(wp), INTENT(inout) :: dheight        ! tree height, [m]
    REAL(wp), INTENT(inout) :: df_lai              ! lai index to estimate tree height




    ! locals
    REAL(wp) :: k_latosa, & ! leaf area to sap wood area ratio [m2/m2]
                sla, & ! specific leaf area [m2/kg C]
                k_rtos, & ! [1/m]
                k_crtos, & ! coarse root to sapwood ratio
                tau_fruit, tau_sapwood, wood_density, frac_sapwood_branch
    REAL(wp) :: leaf2sapwood_mass_ratio, & ! k_latosa / sla / wood_density [m]
                leaf2root_mass_ratio, & 
                leaf2wood_mass_ratio, & ! leaf2sapwood_mass_ratio/ (1 + lctlib_k_crtos)
                frac_alloc_fruit, &
                leaf_carbon_inc, root_carbon_inc, coarse_root_carbon_inc, sap_wood_carbon_inc
    REAL(wp) :: fruit_2_seedbed, sap_2_heartwood
    INTEGER  :: quincy_LctNumber
    ! locals from original calc_Cpools
    REAL(wp) :: NPP_2_greenPool, NPP_2_woodPool, NPP_2_reservePool, NPP_2_rootExudates, Green_2_herbivory
    REAL(wp) :: excess_carbon

    REAL(wp) :: c_green_pot, c_woods_pot, c_reserve_pot        !! values of potential carbon pools i.e. before accounting for N-limitation
    REAL(wp) :: t_tree_mature, t_lai_max, km_lai

    !> init local variables
    leaf_carbon_inc        = 0._wp
    root_carbon_inc        = 0._wp
    coarse_root_carbon_inc = 0._wp
    sap_wood_carbon_inc    = 0._wp
    fruit_2_seedbed        = 0._wp
    sap_2_heartwood        = 0._wp

    t_tree_mature = 0._wp
    t_lai_max     = 0._wp
    km_lai        = 0._wp

    !> retrieve QUINCY lctlib values for turnover time
    CALL convert_jsb_lct_to_quincy(jsb_LctNumber, quincy_LctNumber, &
      &   lctlib_k_latosa = k_latosa, &
      &   lctlib_sla = sla, &
      &   lctlib_wood_density = wood_density, &
      &   lctlib_k_crtos = k_crtos, &
      &   lctlib_tau_fruit = tau_fruit, &
      &   lctlib_tau_sapwood = tau_sapwood, &
      &   lctlib_frac_sapwood_branch = frac_sapwood_branch, &
      &   lctlib_t_tree_mature = t_tree_mature, &
      &   lctlib_t_lai_max = t_lai_max, &
      &   lctlib_km_lai = km_lai
      &   )

    !> retrieve QUINCY 'fixed' allocation coefficient for veg organisms
    SELECT CASE(pheno_type)
    CASE(1, 2, 3)
    ! trees
      leaf2sapwood_mass_ratio = k_latosa / sla / wood_density
      k_rtos = 1.75_wp * (15.4_wp * wood_density) / k_latosa
      leaf2root_mass_ratio    = k_rtos * leaf2sapwood_mass_ratio
      leaf2wood_mass_ratio = leaf2sapwood_mass_ratio / (1._wp + k_crtos)
    CASE(4)
    ! grasses
      leaf2sapwood_mass_ratio = 1._wp / sm2lm_grass
      k_rtos = 0.5_wp / sm2lm_grass
      leaf2root_mass_ratio    = sm2lm_grass * k_rtos
      leaf2wood_mass_ratio = 0._wp ! not applicable
    CASE(5)
    ! crops
      leaf2sapwood_mass_ratio = 1._wp / sm2lm_grass
      k_rtos = 0.5_wp / sm2lm_grass
      leaf2root_mass_ratio    = sm2lm_grass * k_rtos
      leaf2wood_mass_ratio = 0._wp ! not applicable
    END SELECT


    !!
    !! ============== START OF POTENTIAL PLANT-CARBON ALLOCATION =================================
    !!                       (exception: case NPP<0 handles actual allocation
    c_green_pot   = c_green
    c_woods_pot   = c_woods
    c_reserve_pot = c_reserve

    !! Preliminary determination of distribution of NPP to the various pools (may be corrected afterwards according to
    !!    avaialable Nitrogen)

    IF (NPP_pot <= 0.0_wp) THEN                             !! In case of negative or zero NPP

      NPP_2_reservePool =  NPP_pot                          !! .. the plant looses C and
      c_reserve     = c_reserve + NPP_2_reservePool         !! .. we try to take it from the reserve pool
      IF (c_reserve < 0.0_wp) THEN                          !! But if thereby the reserve pool gets negative
        NPP_2_reservePool = NPP_2_reservePool - c_reserve   !! .. we absorb negative NPP only up to zero reserve pool and
        NPP_flux_correction = -c_reserve / sec_per_day      !! .. give the deficit as a flux correction back to the
        c_reserve = 0.0_wp                                  !! .. calling routine. Hence all Carbon from the reserve pool
      END IF                                                !! .. is used up.

      !! For correct handling of nitrogen limitation below (which handles only the case of positive NPP) set NPP to zero:
      NPP_2_reservePool  = 0.0_wp
      NPP_2_greenPool    = 0.0_wp   !! All other transfer rates are zero
      NPP_2_woodPool     = 0.0_wp
      NPP_2_rootExudates = 0.0_wp
      c_reserve_pot = c_reserve    !! remember modified reserve pool
      
      leaf_carbon_inc        = 0._wp
      root_carbon_inc        = 0._wp
      coarse_root_carbon_inc = 0._wp
      sap_wood_carbon_inc    = 0._wp

    ELSE  !! NPP is positive
      ! root exudates and flux to reserve pool is calculated following jsb
      NPP_2_rootExudates= fract_npp_2_exudates * NPP_pot
      NPP_2_reservePool = fract_npp_2_reservePool * NPP_pot
      IF (jsb_allocation) THEN
        ! original JSBACH allocation
        NPP_2_woodPool    = fract_npp_2_woodPool * NPP_pot  !! NPP it is distributed according to predefined
                                                                                   !!    relative fraction
        NPP_2_greenPool   = (1.0_wp - fract_npp_2_woodPool - fract_npp_2_reservePool                 &
             - fract_npp_2_exudates) * NPP_pot

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
        ! update QUINCY pools with fixed allocation distributions
        SELECT CASE(pheno_type)
        CASE(1, 2, 3)
          ! trees
          leaf_carbon_inc = NPP_2_greenPool * 0.45_wp ! 45%
          root_carbon_inc = NPP_2_greenPool * 0.4_wp  ! 40%
          frac_alloc_fruit = 1._wp / 3_wp ! 15%
          coarse_root_carbon_inc = NPP_2_woodPool / 3._wp ! one third
          sap_wood_carbon_inc    = NPP_2_woodPool / 3._wp * 2._wp
        CASE(4, 5)
          ! grasses, crops
          leaf_carbon_inc = NPP_2_greenPool * 0.45_wp ! 45%
          root_carbon_inc = NPP_2_greenPool * 0.4_wp  ! 40%
          frac_alloc_fruit = 1._wp / 3_wp ! 15%
          sap_wood_carbon_inc    = NPP_2_woodPool
        END SELECT


      ELSE
        ! QUINCY allocation
        frac_alloc_fruit = 1._wp / 3_wp ! allocation to fruit is 1/3 of leaf_carbon_inc
        CALL simple_pipe_model(pheno_type, NPP_pot, &
                              c_leaf, c_fineroot, c_sapwood, c_heartwood, c_coarseroot, &
                              leaf_carbon_inc, root_carbon_inc, coarse_root_carbon_inc, sap_wood_carbon_inc, &
                              leaf2sapwood_mass_ratio, k_rtos, leaf2root_mass_ratio, leaf2wood_mass_ratio, k_crtos, frac_alloc_fruit, &
                              tree_height) ! , f_lai, dheight, df_lai


      END IF
    END IF !! NPP > 0 end  

    !! ============== update tree height and f_lai ============================================
    tree_height = tree_height + dheight
    f_lai       = f_lai       + df_lai

    !! ============== use nitrogen from N mobile pool to satisfy (at least partly) for potential carbon fluxes ===
    ! Remarks: Nmobile_2_green + Nmobile_2_wood <= Npool_mobile
    ! Npool_mobile accounts minN_plant_demand by trasfering N into green and wood pool

    !! ============== UPDATE ALL C-POOLS IN JSB & QUINCY =======================================================================
    !> potential C increment of each organism
    
    ! npp to green pools QUINCY
    c_leaf          = c_leaf     + leaf_carbon_inc  - lit_leaf
    c_fineroot      = c_fineroot + root_carbon_inc  - lit_fineroot
    fruit_2_seedbed = MIN(c_fruit * tau_fruit, c_fruit + leaf_carbon_inc * frac_alloc_fruit - lit_fruit)
    c_fruit         = c_fruit    + leaf_carbon_inc * frac_alloc_fruit - lit_fruit - fruit_2_seedbed
    c_seedbed       = c_seedbed  + fruit_2_seedbed  - lit_seedbed
    
    ! npp to wood pools QUINCY
    c_coarseroot   = c_coarseroot + coarse_root_carbon_inc - lit_coarseroot
    sap_2_heartwood = MIN((1._wp - frac_sapwood_branch) * tau_sapwood, c_sapwood + sap_wood_carbon_inc - lit_sapwood)
    c_sapwood      = c_sapwood    + sap_wood_carbon_inc - lit_sapwood - sap_2_heartwood   
    c_heartwood    = c_heartwood  + sap_2_heartwood - lit_heartwood
    
    ! JSBACH pools
    c_green       = c_green   + NPP_2_greenPool
    c_woods       = c_woods   + NPP_2_woodPool
    c_reserve     = c_reserve + NPP_2_reservePool        !!dk removed Nlimit, as no dependent on N


    !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

    IF (NPP_pot > 0.0_wp) THEN
      NPP_act = (NPP_2_rootExudates + NPP_2_greenPool + NPP_2_woodPool + NPP_2_reservePool  &
        & ) / sec_per_day
    ELSE
      NPP_act = NPP_pot / sec_per_day + NPP_flux_correction  ! NPP_pot_yDayMean is negative: NPP_flux_correction is zero
                                                                 ! except c_reserve was < 0.0_wp.
    END IF    

  END SUBROUTINE sb_carbon_allometry

  ! --- convert_lct_to_sla() ---------------------------------------------------------------------------------------------------
  !
  ! find the sla value for JSBACH lctNumber
  !  
  SUBROUTINE convert_jsb_lct_to_quincy(jsb_LctNumber, lctlib_sla, lctlib_cn_leaf, lctlib_np_leaf, &
                 &                     lctlib_tau_leaf, lctlib_tau_fine_root, lctlib_tau_coarse_root, &
                 &                     lctlib_tau_branch, lctlib_tau_wood, lctlib_tau_fruit, lctlib_tau_seedbed, &
                 &                     lctlib_frac_sapwood_branch, lctlib_k_crtos, lctlib_k_latosa, lctlib_wood_density, & 
                 &                     lctlib_t_tree_mature, lctlib_t_lai_max, lctlib_km_lai)
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
    INTEGER, INTENT(in)                       :: jsb_LctNumber              !< lct number from JSBACH lctlib file
    INTEGER, INTENT(out)                      :: quincy_LctNumber           !< corresponding QUINCY lct number
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_sla                 !< corresponding QUINCY lct sla
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_cn_leaf             !< corresponding QUINCY lct cn_leaf
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_np_leaf             !< corresponding QUINCY lct np_leaf
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_k_crtos             !< corresponding QUINCY lct k_crtos
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_k_latosa            !< corresponding QUINCY lct k_latosa
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_leaf            !< corresponding QUINCY lct tau_leaf
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_fine_root       !< corresponding QUINCY lct tau_fine_root
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_coarse_root     !< corresponding QUINCY lct tau_coarse_root
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_branch          !< corresponding QUINCY lct tau_branch
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_wood            !< corresponding QUINCY lct tau_sap_wood
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_fruit           !< corresponding QUINCY lct tau_fruit
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_tau_seedbed         !< corresponding QUINCY lct tau_seed_litter
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_frac_sapwood_branch !< corresponding QUINCY lct frac_sapwood_branch
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_wood_density        !< corresponding QUINCY lct wood_density

    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_t_tree_mature       !< JSBACH tree mature time
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_t_lai_max           !< JSBACH lai saturation time
    REAL(wp), OPTIONAL, INTENT(out)           :: lctlib_km_lai              !< JSBACH km_lai for the height function
    ! ---------------------------
    ! 0.2 Local
    ! --------------------------------------------------------------------
    ! QUINCY LCTLIB lookup arrays
    ! Index: quincy_LctNumber = 1..14
    ! --------------------------------------------------------------------

    REAL(wp), PARAMETER :: q_sla_array(14) = (/ &
        &  9.46_wp,  8.99_wp, 15.40_wp, 15.40_wp,  5.00_wp, 10.09_wp, &
        & 20.12_wp, 19.23_wp, 20.12_wp, 19.23_wp, 20.12_wp, 19.23_wp, &
        & 19.23_wp, 19.23_wp &
        /)

    REAL(wp), PARAMETER :: q_cn_leaf_array(14) = (/ &
        & 16.89_wp, 13.73_wp, 21.32_wp, 21.32_wp, 12.09_wp, 19.37_wp, &
        & 17.84_wp, 14.14_wp, 17.84_wp, 14.14_wp, 17.84_wp, 14.14_wp, &
        & 14.14_wp, 14.14_wp &
        /)

    REAL(wp), PARAMETER :: q_np_leaf_array(14) = (/ &
        & 0.86_wp, 0.84_wp, 1.44_wp, 1.44_wp, 1.23_wp, 1.83_wp, &
        & 1.43_wp, 1.36_wp, 1.43_wp, 1.36_wp, 1.43_wp, 1.36_wp, &
        & 1.36_wp, 1.36_wp &
        /)

    REAL(wp), PARAMETER :: q_k_crtos_array(14) = (/ &
        & 0.1_wp, 0.1_wp, 0.1_wp, 0.1_wp, 0.1_wp, 0.1_wp, &
        & 0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp, &
        & 0.2_wp, 0.2_wp &
        /)

    REAL(wp), PARAMETER :: q_k_latosa_array(14) = (/ &
        & 4000._wp, 4000._wp, 4000._wp, 4000._wp, 4000._wp, 4000._wp, &
        &    0._wp,    0._wp,    0._wp,    0._wp,    0._wp,    0._wp, &
        &    0._wp,    0._wp &
        /)

    REAL(wp), PARAMETER :: q_tau_leaf_array(14) = (/ &
        & 16.83_wp, 15.88_wp,  5.83_wp,  5.83_wp, 39.70_wp,  6.08_wp, &
        &  3.85_wp,  3.85_wp,  3.85_wp,  3.85_wp,  3.85_wp,  3.85_wp, &
        &  3.85_wp,  3.85_wp &
        /) / 12._wp

    REAL(wp), PARAMETER :: q_tau_fine_root_array(14) = (/ &
        & 0.7_wp, 0.7_wp, 0.7_wp, 0.7_wp, 0.7_wp, 0.7_wp, &
        & 0.7_wp, 0.7_wp, 0.7_wp, 0.7_wp, 0.7_wp, 0.7_wp, &
        & 0.7_wp, 0.7_wp &
        /)

    REAL(wp), PARAMETER :: q_tau_coarse_root_array(14) = (/ &
        & 8.0_wp, 8.0_wp, 8.0_wp, 8.0_wp, 8.0_wp, 8.0_wp, &
        & 8.0_wp, 8.0_wp, 8.0_wp, 8.0_wp, 8.0_wp, 8.0_wp, &
        & 8.0_wp, 8.0_wp &
        /)

    REAL(wp), PARAMETER :: q_tau_branch_array(14) = (/ &
        & 10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp, &
        & 10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp, &
        & 10.0_wp, 10.0_wp &
        /)

    REAL(wp), PARAMETER :: q_tau_wood_array(14) = (/ &
        & 40.0_wp, 40.0_wp, 40.0_wp, 40.0_wp, 40.0_wp, 40.0_wp, &
        & 40.0_wp, 40.0_wp, 40.0_wp, 40.0_wp, 40.0_wp, 40.0_wp, &
        & 40.0_wp, 40.0_wp &
        /)

    REAL(wp), PARAMETER :: q_tau_fruit_array(14) = (/ &
        & 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
        & 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
        & 1.0_wp, 1.0_wp &
        /)

    REAL(wp), PARAMETER :: q_tau_seedbed_array(14) = (/ &
        & 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, &
        & 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, &
        & 2.0_wp, 2.0_wp &
        /)

    REAL(wp), PARAMETER :: q_frac_sapwood_branch_array(14) = (/ &
        & 0.05_wp, 0.05_wp, 0.05_wp, 0.05_wp, 0.05_wp, 0.05_wp, &
        & 1.00_wp, 1.00_wp, 1.00_wp, 1.00_wp, 1.00_wp, 1.00_wp, &
        & 1.00_wp, 1.00_wp &
        /)

    REAL(wp), PARAMETER :: q_wood_density_array(14) = (/ &
        & 0.30_wp, 0.30_wp, 0.30_wp, 0.30_wp, 0.30_wp, 0.30_wp, &
        & 0.22_wp, 0.22_wp, 0.22_wp, 0.22_wp, 0.22_wp, 0.22_wp, &
        & 0.22_wp, 0.22_wp &
        /)
    ! --------------------------------------------------------------------
    ! JSBACH LCTLIB lookup arrays
    ! Index: jsb_LctNumber = 1..21
    ! --------------------------------------------------------------------
    REAL(wp), PARAMETER :: j_t_mature_array(21) = (/ &
        & 0._wp, 60._wp, 50._wp, 70._wp, 60._wp, 50._wp, 40._wp, 80._wp, 90._wp, &
        & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
        & 0._wp, 0._wp &
        /)

    REAL(wp), PARAMETER :: j_t_lai_max_array(21) = (/ &
        & 0._wp, 15._wp, 10._wp, 15._wp, 10._wp, 10._wp,  8._wp, 20._wp, 20._wp, &
        & 0._wp,  0._wp,  0._wp,  0._wp,  0._wp,  0._wp,  0._wp,  0._wp,  0._wp, &
        & 0._wp,  0._wp,  0._wp &
        /)

    REAL(wp), PARAMETER :: j_km_lai_array(21) = (/ &
        & 0._wp, 0.6667_wp, 0.6667_wp, 1.0000_wp, 0.6667_wp, 0.6667_wp, 0.6667_wp, &
        & 1.0000_wp, 1.0000_wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
        & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp &
        /)

    !> init out var
    !>
    
    SELECT CASE(jsb_LctNumber)
    CASE(2, 6)
      !> moist broadleaved evergreen
      quincy_LctNumber = 1
    CASE(4)
      !> dry broadleaved evergreen
      quincy_LctNumber = 2
    CASE(3, 5, 10)
      !> rain green broadleaved deciduous
      quincy_LctNumber = 3
    CASE(7, 11)
      !> summer green broadleaved deciduous
      quincy_LctNumber = 4
    CASE(8)
      !> needle-leaved evergreen
      quincy_LctNumber = 5
    CASE(9)
      !> summer green needle-leaved
      quincy_LctNumber = 6
    CASE(12)
      !> C3 grass
      quincy_LctNumber = 7
    CASE(14, 15)
      !> C3 pasture
      quincy_LctNumber = 9
    CASE(19, 20)
      !> C3 crop
      quincy_LctNumber = 11
    CASE(13)
      !> C4 grass
      quincy_LctNumber = 8
    CASE(16)
      !> C4 pasture
      quincy_LctNumber = 10
    CASE(21)
      !> C4 crop
      quincy_LctNumber = 12
    CASE(1, 17, 18)
      !> glacier, Tundra, Swamp
      quincy_LctNumber = 13 
    END SELECT          

    ! --- QUINCY LCTLIB optional values -----------------------------------------
    IF (PRESENT(lctlib_sla))                 lctlib_sla                 = q_sla_array(quincy_LctNumber)
    IF (PRESENT(lctlib_cn_leaf))             lctlib_cn_leaf             = q_cn_leaf_array(quincy_LctNumber)
    IF (PRESENT(lctlib_np_leaf))             lctlib_np_leaf             = q_np_leaf_array(quincy_LctNumber)
    IF (PRESENT(lctlib_k_crtos))             lctlib_k_crtos             = q_k_crtos_array(quincy_LctNumber)
    IF (PRESENT(lctlib_k_latosa))            lctlib_k_latosa            = q_k_latosa_array(quincy_LctNumber)

    IF (PRESENT(lctlib_tau_leaf))            lctlib_tau_leaf            = q_tau_leaf_array(quincy_LctNumber)
    IF (PRESENT(lctlib_tau_fine_root))       lctlib_tau_fine_root       = q_tau_fine_root_array(quincy_LctNumber)
    IF (PRESENT(lctlib_tau_coarse_root))     lctlib_tau_coarse_root     = q_tau_coarse_root_array(quincy_LctNumber)
    IF (PRESENT(lctlib_tau_branch))          lctlib_tau_branch          = q_tau_branch_array(quincy_LctNumber)
    IF (PRESENT(lctlib_tau_wood))            lctlib_tau_wood            = q_tau_wood_array(quincy_LctNumber)
    IF (PRESENT(lctlib_tau_fruit))           lctlib_tau_fruit           = q_tau_fruit_array(quincy_LctNumber)
    IF (PRESENT(lctlib_tau_seedbed))         lctlib_tau_seedbed         = q_tau_seedbed_array(quincy_LctNumber)

    IF (PRESENT(lctlib_frac_sapwood_branch)) lctlib_frac_sapwood_branch = q_frac_sapwood_branch_array(quincy_LctNumber)
    IF (PRESENT(lctlib_wood_density))        lctlib_wood_density        = q_wood_density_array(quincy_LctNumber)
    ! --- JSBACH LCTLIB optional values -----------------------------------------
    IF (PRESENT(lctlib_t_tree_mature))       lctlib_t_tree_mature       = j_t_mature_array(jsb_LctNumber)
    IF (PRESENT(lctlib_t_lai_max))           lctlib_t_lai_max           = j_t_lai_max_array(jsb_LctNumber)
    IF (PRESENT(lctlib_km_lai))              lctlib_km_lai              = j_km_lai_array(jsb_LctNumber)

  END SUBROUTINE convert_jsb_lct_to_quincy
  
  ! --- simple_pipe_model() ---------------------------------------------------------------------------------------------------
  !
  ! Inputs: NPP, L, R, S, C, height, dheight,
  !         r_LR, r_LS, r_LW, frac_alloc_fruit

  ! h_new = height + dheight

  ! R_new = L_new / r_LR
  ! S_new = L_new / r_LS * h_new
  ! C_new = k_crtos * S_new
  ! --> S_new + C_new = L_new / r_LW * h_new
  ! dF = dL * frac_alloc_fruit = (L_new - L) * frac_alloc_fruit
  
  ! NPP + L + R + S + C = L_new + R_new + S_new + C_new + dF
  ! -->                 = L_new + L_new / r_LR + L_new / r_LW * h_new + (L_new - L) * frac_alloc_fruit
  ! NPP + (1+frac_alloc_fruit)*L + R + S + C = L_new * (1 + 1/r_LR + h_new/r_LW + frac_alloc_fruit)

  ! A = 1._wp + 1._wp/r_LR + h_new/r_LW + frac_alloc_fruit
  ! B = NPP + (1+frac_alloc_fruit)*L + R + S + C

  ! L_new = B / A

  ! dL = L_new - L
  ! dR = R_new - R
  ! dS = S_new - S
  ! dC = C_new - C
  ! simple pipe model following quincy fixed allocation scheme  
  SUBROUTINE simple_pipe_model(pheno_type, NPP, &
                              c_leaf, c_fineroot, c_sapwood, c_coarseroot, &
                              d_leaf, d_root, d_croot, d_swood, &
                              r_LS, r_LR, r_LW, k_crtos, frac_alloc_fruit, &
                              h) ! , f_lai, dh, df_lai
    ! USE 

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,  INTENT(in)    :: pheno_type         ! phenology type from JSB lctlib
    REAL(wp), INTENT(in)    :: NPP                ! potential NPP [mol(C)/m^2(canopy)/yr]
    REAL(wp), INTENT(in)    :: c_leaf             ! Leaf carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(in)    :: c_fineroot         ! Fineroot carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(in)    :: c_sapwood          ! SapWood carbon pool (sap and heart wood combined): on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(in)    :: c_coarseroot       ! Coarseroot carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: d_leaf             ! leaf carbon increment [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: d_root             ! fineroot carbon increment [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: d_croot            ! Coarseroot carbon increment [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: d_wood            ! sapwood carbon increment [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(in)    :: r_LS               ! leaf2sapwood_mass_ratio, [m]
    REAL(wp), INTENT(in)    :: r_LR               ! leaf2root_mass_ratio, [-]
    REAL(wp), INTENT(in)    :: r_LW               ! leaf2wood_mass_ratio, [m]
    REAL(wp), INTENT(in)    :: frac_alloc_fruit   ! ratio between d_fruit and d_leaf
    REAL(wp), INTENT(in)    :: h                  ! tree height, [m]
    ! REAL(wp), INTENT(in)    :: f_lai              ! lai index to estimate tree height
    ! REAL(wp), INTENT(inout) :: dh                 ! tree height increment, [m]
    ! REAL(wp), INTENT(in)    :: df_lai             ! lai index increment

    ! to_do: 
    ! df_lai and dh should be updated in update_phenology_logrop SUBROUTINE and passed to SB_CARBON as inputs

    ! locals
    REAL(wp) :: c_leaf_new, c_fineroot_new, c_sapwood_new, c_coarseroot_new, c_wood_new
    REAL(wp) :: numerator, denominator
    
    ! init locals
    d_leaf  = 0._wp
    d_root  = 0._wp
    d_croot = 0._wp
    d_wood  = 0._wp
    
    
    SELECT CASE(pheno_type)
    !> trees
    !!
    CASE (1, 2, 3)
      numerator   = NPP + (1._wp+frac_alloc_fruit)*c_leaf + c_fineroot + c_sapwood + c_coarseroot
      denominator = 1._wp + 1._wp/r_LR + h/r_LW + frac_alloc_fruit
      ! leaf
      c_leaf_new  = numerator / denominator
      d_leaf      = c_leaf_new - c_leaf
      ! SapWood
      c_sapwood_new = c_leaf_new / r_LS * h
      d_wood        = c_sapwood_new - c_sapwood
      IF (d_wood > 0._wp) Then
        d_croot = MAX(c_sapwood_new * k_crtos - c_coarseroot, 0.0_wp)
        d_root  = MAX(c_leaf_new / r_LR - c_root, 0.0_wp)
      ELSE! special case of no wood growth (phenology etc.): redistribute available carbon
        numerator   = NPP + (1._wp+frac_alloc_fruit)*c_leaf + c_fineroot
        denominator = 1._wp + 1._wp/r_LR + frac_alloc_fruit
        c_leaf_new  = numerator / denominator
        d_leaf      = MAX(c_leaf_new - c_leaf, 0.0_wp)
        d_root      = MAX(0.0_wp, MIN(NPP - (1._wp+frac_alloc_fruit) * d_leaf, (d_leaf + c_leaf) / r_LR - c_root))
        d_wood      = 0.0_wp
        d_croot     = 0.0_wp
      END IF

    !> grass & crops
    !!
    CASE(4, 5)
      numerator   = NPP + (1+frac_alloc_fruit)*c_leaf + c_fineroot + c_sapwood
      denominator = 1._wp + 1._wp/r_LR + h/r_LS + frac_alloc_fruit
      c_leaf_new  = numerator / denominator
      d_leaf      = MAX(c_leaf_new - c_leaf, 0.0_wp)

      d_root      = MAX(MIN(NPP-(1._wp+frac_alloc_fruit) * d_leaf,
                            (c_leaf + d_leaf)/r_LR - c_root),0.0_wp)
      d_wood      = (d_leaf + c_leaf) / r_LS - c_sapwood
      d_croot     = 0.0_wp
    CASE DEFAULT
      !! Write statement
      !!
      WRITE(message_text,'(a,i4)') 'invalid pheno type ',pheno_type
      CALL finish("simple_pipe_model",message_text)
    END SELECT

  END SUBROUTINE simple_pipe_model
  
  ! --- calc_sb_treeheight() ---------------------------------------------------------------------------------------------------
  !
  ! calculate tree height for sb_carbon simple pipe model
  ! 
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_sb_treeheight(jsb_LctNumber, dtime, LAI_old, LAI, MaxLAI, MaxHeight f_lai, df_lai, h, dh)

    USE mo_carbon_constants,    ONLY: sec_per_year
    USE mo_jsb_math_constants,  ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER, INTENT(in)     :: jsb_LctNumber      ! LctNumber from jsb lctlib
    REAL(wp), INTENT(in)    :: dtime              ! dtime
    REAL(wp), INTENT(in)    :: LAI_old            ! LAI from previous dt
    REAL(wp), INTENT(in)    :: LAI                ! LAI after updating
    REAL(wp), INTENT(in)    :: MaxLAI             ! lctlib maximim LAI
    REAL(wp), INTENT(in)    :: MaxHeight          ! lctlib maximim height
    REAL(wp), INTENT(inout) :: h                  ! tree height, [m]
    REAL(wp), INTENT(inout) :: f_lai              ! lai index to estimate tree height
    REAL(wp), INTENT(inout) :: dh                 ! tree height increment, [m]
    REAL(wp), INTENT(inout) :: df_lai             ! lai index increment

    ! locals
    REAL(wp) :: f_lai_realtime, & ! realtime f_lai calculated by LAI/MaxLAI
                t_lai_inc         ! years to reach f_lai phase, canopy closure or maximum height
    REAL(wp) :: t_tree_mature, t_lai_max, km_lai, h_new

    !> init local variables
    f_lai_realtime = 0._wp
    t_lai_inc      = 0._wp

    t_tree_mature = 0._wp
    t_lai_max     = 0._wp
    km_lai        = 0._wp
    !> init output
    df_lai = 0._wp
    dh     = 0._wp
    h_new  = h

    !> retrieve QUINCY lctlib values for turnover time
    CALL convert_jsb_lct_to_quincy(jsb_LctNumber, quincy_LctNumber, &
      &   lctlib_t_tree_mature = t_tree_mature, &
      &   lctlib_t_lai_max = t_lai_max, &
      &   lctlib_km_lai = km_lai
      &   ) 
    
    f_lai_realtime = LAI / MaxLAI

    !> calculate teh change of f_lai based on canopy closure status
    !>
    IF ( (f_lai - km_lai) > eps8 ) THEN
      ! canopy closure, slow growing
      t_lai_inc = t_tree_mature - t_lai_max
      IF ( (LAI - LAI_old) > eps8 ) THEN
        !> ! growing phase, height increase
        df_lai = MaxLAI/t_lai_inc/sec_per_year*dtime
      END IF
    ELSE
      ! fast growing early phase
      t_lai_inc = t_lai_max
      IF ( (LAI - LAI_old) > eps8) THEN
        !> ! growing phase, height increase
        IF ( f_lai_realtime > f_lai) THEN
          df_lai = f_lai_realtime - f_lai
        ELSE
          df_lai = MaxLAI/t_lai_inc/sec_per_year*dtime
        END IF
      END IF      
    END IF

    !> calculate new tree height
    IF (df_lai > 0._wp) THEN
      h_new = MaxHeight * (f_lai + df_lai)/(km_lai + f_lai + df_lai)
      dh    = h_new - h
    END if
    
    ! Update h and f_lai
    f_lai = f_lai + df_lai
    h     = h + dh
    END IF

  
  END SUBROUTINE calc_sb_treeheight
  
END MODULE mo_sb_carbon_process
