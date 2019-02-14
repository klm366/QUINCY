!>
!! @brief Contains the memory class for the carbon process.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-6)
!!
MODULE mo_carbon_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  !USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  !USE mo_jsb_lct_class,          ONLY: VEG_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d !, t_jsb_var_real3d

  ! Use of prcesses in this module
  dsl4jsb_Use_processes PHENO_

  ! Use process configurations
  dsl4jsb_Use_config(PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_carbon_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 150    !PHS changed from 110 to 120

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_carbon_memory
    ! 1d reals:
    !REAL(wp) ::            &
    TYPE(t_jsb_var_real2d) ::     &
       N_pseudo_temp,       & ! Constants used for updating the pseudo-15-day mean temperature
       F_pseudo_temp,       & ! Constants used for updating the pseudo-15-day mean temperature
       N_pseudo_precip,     & ! Constants used for updating the pseudo-15-day mean precipitation
       F_pseudo_precip        ! Constants used for updating the pseudo-15-day mean precipitation

    ! cbalance_type variables
    TYPE(t_jsb_var_real2d) ::     &
      pseudo_temp,                & ! a proxy for 15-day mean air temperature, [K] needed with yasso
      pseudo_temp_yDay,           & !
      pseudo_precip,              & ! a proxy for 15-day mean preciptation rate [Kg/m^2 s] needed with yasso
      pseudo_precip_yDay,         & !
      topSoilTemp,                & ! R: in the moment it is still better to leaf topSoilTemp in the CARBON process as it is
                                    !    only used here.
      alpha                         ! R: in the moment it is better to leaf alpha in the CARBON process as it is only used here
                                    !    However, as we used the term alpha also in mo_assimi_constants.f90 for the efficiency of
                                    !    photo capture, we should also use another term.


    TYPE(t_jsb_var_real2d) ::     &
      cconservation_calcCpools,   &
      cflux_c_greenwood_2_litter, & ! R: From JS3: cbalance_diag_type. Carbon flux from the vegetation to the litter pools
                                    !    [mol(C)/m^2(canopy) s]
      c_green,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      c_woods,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      c_reserve,              & ! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive
                                    ! .. bad times[mol(C)/m^2(canopy)]
!@todo why is the c_crop_harvest realised as a density and not as a per tile variable in the first place?
! ++ only required for crop pfts, or?
      c_crop_harvest,         & ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
!@todo: the only variable of these which seems to be still required is c_litter_green_ag - rest could be local/removed
!       (or at least commented out as recently done for c_slow) (including ta variables!)
      c_litter_green_ag,      & ! C-pool for above ground non-woody litter especially from leaves
                                    !     [mol(C)/m^2(canopy)]
      c_litter_green_bg,      & ! C-pool for below ground non-woody litter especially from fine roots
                                    !     [mol(C)/m^2(canopy)]
      c_litter_wood_ag,       & ! C-pool for litter originating from above ground woody parts of the plants
                                    !     [mol(C)/m^2(canopy)]
      c_litter_wood_bg,       & ! C-pool for litter originating from below ground woody parts of the plants
                                    !     [mol(C)/m^2(canopy)]
     !CBAL c_slow,                 & ! C-pool for below ground organic material (coming from the litter pools)
     !CBAL                              !     [mol(C)/m^2(canopy)]
     !SIZE CLASS 1
      c_acid_ag1,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_ag1,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_ag1,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_ag1,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_acid_bg1,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_bg1,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_bg1,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_bg1,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_humus_1,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(canopy)]
     !SIZE CLASS 2
      c_acid_ag2,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_ag2,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_ag2,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_ag2,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_acid_bg2,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      c_water_bg2,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      c_ethanol_bg2,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      c_nonsoluble_bg2,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      c_humus_2,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(canopy)]
      c_bg_sum,              & ! Sum of belowground Yasso pools  [mol(C)^2/m(canopy)]
      c_ag_sum_1,            & ! Sum of aboveground green Yasso pools  [mol(C)^2/m(canopy)]

!@todo: why? its not used (i.e. not even caluclated, only the ta variable is calculated)
!      C_sum_natural,              & ! R: added
      !
      LAI_sum,                    & ! used to accumulate LAI over a day
      NPP_pot_sum,                    & ! used to accumulated NPP-Rate over a day
      GPP_sum,                    & ! used to accumulated GPP-Rate over a day
      topSoilTemp_sum,            & ! used to accumulated topSoilTemp over a day. R: Nicht notwendig, kommt nur in Output.
      alpha_sum,                  & ! used to accumulated water stress coefficient alpha over a day
      !runoff_sum,                 & ! used to accumulate runoff over a day.  R: Nicht notwendig in JS4.
      soil_respiration,           & ! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(canopy)].
                                    ! Without N limitation!
      NPP_flux_correction,        & ! Daily updated flux correction from yesterdays carbon balance
                                    ! [mol(CO2)/m^2(canopy) s]
                                    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                    ! correction arises either because otherwise the reserve pool would get
                                    ! negative (positive correction).
      excess_NPP,                 & ! That part of NPP that because of structural limits could not be allocated in
                                    ! carbon pools [mol(CO2)/m^2(canopy) s]
      LAI_yyDayMean,              & ! previous previous day mean
      LAI_yDayMean,               & ! previous day mean
      NPP_pot_yDayMean,           & ! mean value of NPP-Rate yesterday (from NPP_pot_sum()) [mol(CO2)/(m^2(canopy) s)]
      NPP_act_yDayMean,           & ! mean value of actual NPP-Rate after N-limitation and excess carbon drop.
                                    ! "yDayMean" because the actual NPP depends on/is calculated from results of the
                                    ! previous day.
      GPP_yDayMean,               & !
      topSoilTemp_yDayMean,       & !
      alpha_yDayMean,             & !
      root_exudates,              & ! Total root exudates entering to the litter green pools [mol(C)/m^2(canopy) s]
      cflux_c_green_2_herb,            & ! Total carbon flux from C pool green to herbivory.
      cflux_herb_2_littergreen,         & ! ..part of cflux_c_green_2_herb that goes into green litter
      cflux_herb_2_atm,      & ! ..part of cflux_c_green_2_herb that goes into atmosphere
      !
      cflux_dist_greenreserve_2_soil,  & ! R: should be later renamed to C_2_GreenLitterPools
      cflux_dist_woods_2_soil,    & ! R: should be later renamed to C_2_WoodLitterPools
      cflux_fire_all_2_atm,       & ! R: should be later renamed to C_2_atmos
                                    ! Carbon immedieatly released by fire per tile [mol(C)m-2(subtile)]
      ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
      !fract_litter_wood_new,       & ! new fraction in above ground wood litter pool

  !JN: variables required when l_forestRegrowth = true, i.e. when maxLAI is derived from biomass via allometric relationships
      current_max_green,          & ! Running maximum of the green pool [mol(C) m-2(canopy)]
                                    ! since last maxLAI calculation
      veg_carbon_at_max_green       ! Total vegetation carbon [mol(C) m-2(canopy)]
                                    ! at the time of the green pool maximum (current_max_green)

   ! Output on tile areas:
   TYPE(t_jsb_var_real2d) ::         &
      cflux_c_greenwood_2_litter_ta, & ! R: From JS3: cbalance_diag_type. Carbon flux from the vegetation to the litter pools
                                       !    [mol(C)/m^2(tile area) s]
      c_green_ta,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts
                                       ! .. of vegetation [mol(C)/m^2(tile area)]
      c_woods_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      c_reserve_ta,              & ! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive
                                       ! .. bad times[mol(C)/m^2(tile area)]
      c_crop_harvest_ta,         & ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]
      c_litter_green_ag_ta,      & ! C-pool for above ground non-woody litter especially from leaves
                                       !     [mol(C)/m^2(tile area)]
      c_litter_green_bg_ta,      & ! C-pool for below ground non-woody litter especially from fine roots
                                       !     [mol(C)/m^2(tile area)]
      c_litter_wood_ag_ta,       & ! C-pool for litter originating from above ground woody parts of the plants
                                       !     [mol(C)/m^2(tile area)]
      c_litter_wood_bg_ta,       & ! C-pool for litter originating from below ground woody parts of the plants
                                       !     [mol(C)/m^2(tile area)]
      !CBAL c_slow_ta,                 & ! C-pool for below ground organic material (coming from the litter pools)
      !CBAL                              !     [mol(C)/m^2(tile area)]
     !SIZE CLASS 1
      c_acid_ag1_ta,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_ag1_ta,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_ag1_ta,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_ag1_ta,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_acid_bg1_ta,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_bg1_ta,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_bg1_ta,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_bg1_ta,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_humus_1_ta,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(tile area)]
     !SIZE CLASS 2
      c_acid_ag2_ta,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_ag2_ta,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_ag2_ta,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_ag2_ta,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_acid_bg2_ta,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      c_water_bg2_ta,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      c_ethanol_bg2_ta,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      c_nonsoluble_bg2_ta,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      c_humus_2_ta,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(tile area)]
!@todo: why? its not used!
!      c_bg_sum_ta,              & ! Sum of belowground Yasso pools  [mol(C)^2/m(tile area)]
      C_sum_natural_ta,              & ! R: added
      !
      !@todo: NPP_pot_sum_ta, GPP_sum_ta, LAI_yyDayMean_ta and LAI_yDayMean_ta are never calculated
!      NPP_pot_sum_ta,                    & ! used to accumulated NPP-Rate over a day
!      GPP_sum_ta,                    & ! used to accumulated GPP-Rate over a day
      soil_respiration_ta,           & ! mean daily rate of heterotrophic (soil) respiration [mol(CO2)/m^2(tile area)].
                                       ! Without N limitation!
      NPP_flux_correction_ta,        & ! Daily updated flux correction from yesterdays carbon balance
                                    ! [mol(CO2)/m^2(tile area) s]
                                    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                    ! correction arises either because otherwise the reserve pool would get
                                    ! negative (positive correction).
      excess_NPP_ta,                 & ! That part of NPP that because of structural limits could not be allocated in
                                    ! carbon pools [mol(CO2)/m^2(tile area) s]
!      LAI_yyDayMean_ta,              & ! previous previous day mean
!      LAI_yDayMean_ta,               & ! previous day mean
      NPP_pot_yDayMean_ta,           & ! mean value of NPP-Rate yesterday (from NPP_pot_sum()) [mol(CO2)/(m^2(tile area) s)]
      NPP_act_yDayMean_ta,           & ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation.
                                       ! Actual NPP (after N-limitation and) excess carbon drop.
                                       ! [mol(CO2)/(m^2(tile area) s)]
      GPP_yDayMean_ta,               & !
      root_exudates_ta,              & ! Total root exudates entering to the litter green pools [mol(C)/m^2(tile area) s]
      cflux_c_green_2_herb_ta,            & ! Total carbon flux from C pool green to herbivory.
      !cflux_herb_2_littergreen_ta,         & ! ..part of cflux_c_green_2_herb that goes into green litter
      !cflux_herb_2_atm_ta,      & ! ..part of cflux_c_green_2_herb that goes into atmosphere
      cflux_dist_greenreserve_2_soil_ta,  & ! R: should be later renamed to C_2_GreenLitterPools_ta
      cflux_dist_woods_2_soil_ta,   & ! R: should be later renamed to C_2_WoodLitterPools_ta
      cflux_fire_all_2_atm_ta,             & ! R: should be later renamed to C_2_atmos_ta
                                    ! Carbon immedieatly released by fire per tile [mol(C)m-2(subtile)]

      ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
      fract_litter_wood_new_ta,      & ! new fraction in above ground wood litter pool
      ! former (JS3) local variables. R: nur auf tile area benötigt!
      cflux_npp_2_atm_ta,            & ! [kg(CO2)/m^2 s]  tile averages of net CO2 fluxes between
      cflux_soilresp_2_atm_ta,       & !  .. biosphere (due to NPP, soil respiration and
      cflux_herb_2_atm_ta,      & !  .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)]
                                       ! .. (positive for emission to atmosphere,
                                       ! ..  negative for absorption by biosphere)
      cflux_npp_2_atm_yday_ta          ! day CO2 flux from actual NPP, required for C conservation test
      
   ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
   ! Variavles for JeDi carbon allocation 
   
   TYPE(t_jsb_var_real2d) ::         &
       c_storage,                    &
       c_storage_temp,               &
       c_seed,                       &
       temp_c_seed,                  &
       c_stem,                       &
       c_leaves,                     &
       c_fnroot,                     &
       c_csroot,                      &
        c_total                       
   TYPE(t_jsb_var_real2d) ::            &
       c_storage_ta,                    &
       c_seed_ta,                       &
       c_stem_ta,                       &
       c_leaves_ta,                     &
       c_fnroot_ta,                     &
       c_csroot_ta,                    &
       c_total_ta
       
    TYPE(t_jsb_var_real2d) ::         &
        ff_grow_T,                     &
        ff_grow_G,                     &
        ff_grow_W,                     &
        ff_npp,                        &
        Flg_die,                       &
        Flg_senesc,                    &
        alpha_jedi,                    &
        MaintRes,                      &
        MaintRes_sum,                  &
        MaintRes_yDayMean,             &
        GrowthResp,                    &
        AutoRes ,                      &
        C_balan,                       &
        C_delta,                       &
        flx_germ
    
   ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  


  CONTAINS
    PROCEDURE :: Init => Init_carbon_memory
  END TYPE t_carbon_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_memory_class'

CONTAINS

  SUBROUTINE Init_carbon_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: jsb_add_var, BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_carbon_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    dsl4jsb_Def_config(PHENO_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface !, canopy_layer        ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_carbon_memory'
    IF (lct_ids(1) > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    !canopy_layer => Get_vgrid('canopy_layer')

    dsl4jsb_Get_config(PHENO_)

    IF (.NOT. ASSOCIATED(mem%pseudo_temp%ptr)) THEN

! R: only if we would want to have these variables in 2d (see mo_carbon_init.f90 for their initialization):
      CALL jsb_add_var( mem, 'n_pseudo_temp', mem%N_pseudo_temp, &
        & hgrid, surface,                                                                              &
        & t_cf('N_pseudo_temp', 'K', 'Normalization for computing pseudo-15-day mean air temperature'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'f_pseudo_temp', mem%F_pseudo_temp, &
        & hgrid, surface,                                                                              &
        & t_cf('F_pseudo_temp', '-', 'Exponential factor used for updating the pseudo-15-day mean temperature'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'n_pseudo_precip', mem%N_pseudo_precip, &
        & hgrid, surface,                                                                              &
        & t_cf('N_pseudo_precip', 'Kg m-2 s-1', 'Normalization for computing pseudo-15-day mean precipitation rate'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL jsb_add_var( mem, 'f_pseudo_precip', mem%F_pseudo_precip, &
        & hgrid, surface,                                                                              &
        & t_cf('F_pseudo_precip', '-', 'Exponential factor used for updating the pseudo-15-day mean precip'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp )


      CALL jsb_add_var( mem, 'pseudo_temp', mem%pseudo_temp, &
        & hgrid, surface,                                                                              &
        & t_cf('pseudo_temp', 'K', 'Pseudo-mean air temperature over previous 15 days'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=MEDIUM,                                                                         &
        & lrestart=.TRUE., initval_r=285.0_wp ) ! initval taken from JSBACH3, arbitary.

      ! R: In the output it should be the same as pseudo_temp, I do not understand why it was written out in JS3.
      CALL jsb_add_var( mem, 'pseudo_temp_yDay', mem%pseudo_temp_yDay, &
        & hgrid, surface,                                                                                &
        & t_cf('pseudo_temp_yDay', 'K', 'Yesterdays pseudo-mean air temperature over previous 15 days'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
        & prefix, suffix,                                                                                &
        & lrestart=.TRUE., initval_r=285.0_wp ) ! initval taken from JSBACH3, arbitary.

      CALL jsb_add_var( mem, 'pseudo_precip', mem%pseudo_precip, &
        & hgrid, surface,                                                                              &
        & t_cf('pseudo_precip', 'Kg m-2 s-1', 'Pseudo-mean precipitation rate over previous 15 days'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=MEDIUM,                                                                         &
        & lrestart=.TRUE., initval_r=3.0E-5_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3, arbitary.

      ! R: In the output it should be the same as pseudo_precip, I do not understand why it was written out in JS3.
      CALL jsb_add_var( mem, 'pseudo_precip_yDay', mem%pseudo_precip_yDay,                             &
        & hgrid, surface,                                                                              &
        & t_cf('pseudo_precip_yDay', 'Kg m-2 s-1', 'Yesterdays pseudo-mean precipitation rate over previous 15 days'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & lrestart=.TRUE., initval_r=3.0E-5_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3, arbitary.

!@todo: relict from CBALANCE, (currently) not required for YASSO
      CALL jsb_add_var( mem, 'topsoiltemp', mem%topSoilTemp,                                           &
        & hgrid, surface,                                                                              &
        & t_cf('topSoilTemp', 'K', 'Top soil temperature.'),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=FULL,                                                                           &
        & lrestart=.TRUE., initval_r=283.15_wp ) ! I "invented" this initval.

!@todo: relict from CBALANCE, (currently) not required for YASSO
      CALL jsb_add_var( mem, 'alpha', mem%alpha,                                                       &
        & hgrid, surface,                                                                              &
        & t_cf('alpha', '-', 'Water Stress Coefficient. Relative soil moisture.'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=FULL,                                                                           &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3's alpha_yDayMean

!@todo: relict from CBALANCE?
      ! Carbon variables on canopy area
      CALL jsb_add_var( mem, 'c_greenwood2litter', mem%cflux_c_greenwood_2_litter,                       &
        & hgrid, surface,                                                                                &
        & t_cf('c_greenwood2litter', 'mol(CO2) m-2(canopy) s-1', 'Total litter flux entering the soil pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
        & prefix, suffix,                                                                                &
        & output_level=FULL,                                                                             &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_green', mem%c_green,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_green', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=FULL,                                                                             &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_woods', mem%c_woods,                                                &
        & hgrid, surface,                                                                                   &
        & t_cf('c_woods', 'mol(C) m-2(canopy)', 'C-Pool for structural material/parts of vegetation.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
        & prefix, suffix,                                                                                   &
        & output_level=FULL,                                                                             &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_reserve', mem%c_reserve,                                                             &
        & hgrid, surface,                                                                                                    &
        & t_cf('c_reserve', 'mol(C) m-2(canopy)', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                                 &
        & prefix, suffix,                                                                                                    &
        & output_level=FULL,                                                                             &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_crop_harvest', mem%c_crop_harvest,                               &
        & hgrid, surface,                                                                                &
        & t_cf('c_crop_harvest', 'mol(C) m-2(canopy)', 'C-Pool for Harvested Material from Crops.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
        & prefix, suffix,                                                                                &
        & output_level=FULL,                                                                             &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_litter_green_ag', mem%c_litter_green_ag,                                     &
        & hgrid, surface,                                                                                            &
        & t_cf('c_litter_green_ag', 'mol(C) m-2(canopy)', 'C-Pool for above ground non-woody litter (leaves).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval in JSBACH3 is 60.0_wp

      CALL jsb_add_var( mem, 'c_litter_green_bg', mem%c_litter_green_bg,                                         &
        & hgrid, surface,                                                                                                &
        & t_cf('c_litter_green_bg', 'mol(C) m-2(canopy)', 'C-Pool for below ground non-woody litter (fine roots).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & output_level=FULL,                                                                   &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval in JSBACH3 is 60.0_wp

      CALL jsb_add_var( mem, 'c_litter_wood_ag', mem%c_litter_wood_ag,                         &
        & hgrid, surface,                                                                              &
        & t_cf('c_litter_wood_ag', 'mol(C) m-2(canopy)', 'C-Pool for above ground woody litter.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=FULL,                                                  &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval in JSBACH3 is 180.0_wp

      CALL jsb_add_var( mem, 'c_litter_wood_bg', mem%c_litter_wood_bg,                         &
        & hgrid, surface,                                                                              &
        & t_cf('c_litter_wood_bg', 'mol(C) m-2(canopy)', 'C-Pool for below ground woody litter.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=FULL,                        &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval in JSBACH3 is 180.0_wp

!CBAL       CALL jsb_add_var( mem, 'c_slow', mem%c_slow,                                                 &
!CBAL         & hgrid, surface,                                                                                  &
!CBAL         & t_cf('c_slow', 'mol(C) m-2(canopy)', 'C-Pool for slowly respirated soil organic material.'), &
!CBAL         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                               &
!CBAL         & prefix, suffix,                                                                                  &
!CBAL         & output_level=FULL,                                                                             &
!CBAL         & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3 from cbalance%c_slow.
!CBAL                                               ! Note for cbalone this value would be:
!CBAL                                               ! mo_cbalone_memory.f90:247: cbalance%c_slow(:,:) = 0.0_dp

      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      ! Variavles for JeDi carbon allocation 
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
      CALL jsb_add_var( mem, 'c_storage', mem%c_storage,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_storage', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.FALSE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
      
      CALL jsb_add_var( mem, 'c_storage_temp', mem%c_storage_temp,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_storage_temp', 'mol(C) m-2(canopy)', 'Storage pool (before allocation )for checking.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
        
      CALL jsb_add_var( mem, 'c_seed', mem%c_seed,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_seed', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.FALSE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )
        
      CALL jsb_add_var( mem, 'temp_c_seed', mem%temp_c_seed,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('temp_c_seed', 'mol(C) m-2(canopy)', 'C-Pool of seed.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )  

      CALL jsb_add_var( mem, 'c_stem', mem%c_stem,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_stem', 'mol(C) m-2(canopy)', 'C-Pool of stem.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.FALSE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL jsb_add_var( mem, 'c_csroot', mem%c_csroot,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_csroot', 'mol(C) m-2(canopy)', 'C-Pool of coarse root.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.FALSE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL jsb_add_var( mem, 'c_fnroot', mem%c_fnroot,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_fnroot', 'mol(C) m-2(canopy)', 'C-Pool of fine root.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.FALSE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL jsb_add_var( mem, 'c_leaves', mem%c_leaves,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_leaves', 'mol(C) m-2(canopy)', 'C-Pool of leaves'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.FALSE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL jsb_add_var( mem, 'c_total', mem%c_total,                              &
        & hgrid, surface,                                                                        &
        & t_cf('c_total','mol(C) m-2(canopy)' ,'Total carbon of six pools'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   

      CALL jsb_add_var( mem, 'c_storage_ta', mem%c_storage_ta,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_storage_ta', 'mol(C) m-2', 'C-Pool for green parts of vegetation.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
        
      CALL jsb_add_var( mem, 'c_seed_ta', mem%c_seed_ta,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_seed_ta', 'mol(C) m-2', 'C-Pool of seed.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL jsb_add_var( mem, 'c_stem_ta', mem%c_stem_ta,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_stem_ta', 'mol(C) m-2', 'C-Pool of stem.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL jsb_add_var( mem, 'c_csroot_ta', mem%c_csroot_ta,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_csroot_ta', 'mol(C) m-2' ,'C-Pool of coarse root.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

      CALL jsb_add_var( mem, 'c_fnroot_ta', mem%c_fnroot_ta,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_fnroot_ta','mol(C) m-2', 'C-Pool of fine root.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
        
      CALL jsb_add_var( mem, 'c_leaves_ta', mem%c_leaves_ta,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('c_leaves_ta', 'mol(C) m-2', 'C-Pool od leaves.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
      
      CALL jsb_add_var( mem, 'c_total_ta', mem%c_total_ta,                              &
        & hgrid, surface,                                                                        &
        & t_cf('c_total_ta','mol(C) m-2','Total carbon of six pools'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
        
      CALL jsb_add_var( mem, 'ff_grow_T', mem%ff_grow_T,                                      &
        & hgrid, surface,                                                                        &
        & t_cf('ff_grow_T', 'Growing condition respect to air temperature.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. ) 
        
      CALL jsb_add_var( mem, 'ff_grow_G', mem%ff_grow_G,                                      &
        & hgrid, surface,                                                                        &
        & t_cf('ff_grow_G', 'Germination condition'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
       
      CALL jsb_add_var( mem, 'ff_grow_W', mem%ff_grow_W,                                      &
        & hgrid, surface,                                                                        &
        & t_cf('ff_grow_W', 'Growing condition respect to soil wetness.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
      
      CALL jsb_add_var( mem, 'ff_npp', mem%ff_npp,                                              &
        & hgrid, surface,                                                                        &
        & t_cf('ff_npp', 'Senescence condition respect to net primary production.'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
      
      CALL jsb_add_var( mem, 'flx_germ', mem%flx_germ,                              &
        & hgrid, surface,                                                                        &
        & t_cf('flx_germ', 'carbon flux from germination '),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
      
      CALL jsb_add_var( mem, 'Flg_die', mem%Flg_die,                                              &
        & hgrid, surface,                                                                        &
        & t_cf('Flg_die', 'Dead'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
    
      CALL jsb_add_var( mem, 'Flg_senesc', mem%Flg_senesc,                                              &
        & hgrid, surface,                                                                        &
        & t_cf('Flg_senesc', 'Senescence flag'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
        
      CALL jsb_add_var( mem, 'alpha_jedi', mem%alpha_jedi,                              &
        & hgrid, surface,                                                                        &
        & t_cf('alpha_jedi', 'relalive soil wetness for JeDi'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
          
      
      CALL jsb_add_var( mem, 'MaintRes', mem%MaintRes,                              &
        & hgrid, surface,                                                                        &
        & t_cf('MaintRes', 'Maintenance Respiration'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.False.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
        
      CALL jsb_add_var( mem, 'MaintRes_sum', mem%MaintRes_sum,                              &
        & hgrid, surface,                                                                        &
        & t_cf('MaintRes_sum', 'Mean maintenance Respiration of previous day'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.False.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )      
        
      CALL jsb_add_var( mem, 'MaintRes_ydaymean', mem%MaintRes_yDayMean,                              &
        & hgrid, surface,                                                                        &
        & t_cf('MaintRes_ydaymean', 'Mean maintenance Respiration of previous day'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
      
      CALL jsb_add_var( mem, 'GrowthResp', mem%GrowthResp,                              &
        & hgrid, surface,                                                                        &
        & t_cf('GrowthResp', 'Growth respiration'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
      
      
      CALL jsb_add_var( mem, 'AutoRes', mem%AutoRes,                              &
        & hgrid, surface,                                                                        &
        & t_cf('AutoRes', 'Autotrophic respiration'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
       
      CALL jsb_add_var( mem, 'C_balan', mem%C_balan,                              &
        & hgrid, surface,                                                                        &
        & t_cf('C_balan', 'Carbon balance check'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
        
       CALL jsb_add_var( mem, 'C_delta', mem%C_delta,                              &
        & hgrid, surface,                                                                        &
        & t_cf('C_delta', 'Carbon delta check'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=BASIC,loutput=.TRUE.,                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
        
        
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

      ! Yasso 1
      CALL jsb_add_var( mem, 'c_acid_ag1', mem%c_acid_ag1,                                               &
        & hgrid, surface,                                                                                          &
        & t_cf('c_acid_ag1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_ag1', mem%c_water_ag1,                                               &
        & hgrid, surface,                                                                                            &
        & t_cf('c_water_ag1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_ag1', mem%c_ethanol_ag1,                                               &
        & hgrid, surface,                                                                                                &
        & t_cf('c_ethanol_ag1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_ag1', mem%c_nonsoluble_ag1,                                        &
        & hgrid, surface,                                                                                               &
        & t_cf('c_nonsoluble_ag1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_acid_bg1', mem%c_acid_bg1,                                               &
        & hgrid, surface,                                                                                          &
        & t_cf('c_acid_bg1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=FULL,                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_bg1', mem%c_water_bg1,                                               &
        & hgrid, surface,                                                                                            &
        & t_cf('c_water_bg1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_bg1', mem%c_ethanol_bg1,                                              &
        & hgrid, surface,                                                                                               &
        & t_cf('c_ethanol_bg1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_bg1', mem%c_nonsoluble_bg1,                                       &
        & hgrid, surface,                                                                                              &
        & t_cf('c_nonsoluble_bg1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
        & prefix, suffix,                                                                                              &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_humus_1', mem%c_humus_1,                     &
        & hgrid, surface,                                                              &
        & t_cf('c_humus_1', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      ! Yasso 2
      CALL jsb_add_var( mem, 'c_acid_ag2', mem%c_acid_ag2,                                               &
        & hgrid, surface,                                                                                          &
        & t_cf('c_acid_ag2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_ag2', mem%c_water_ag2,                                               &
        & hgrid, surface,                                                                                            &
        & t_cf('c_water_ag2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_ag2', mem%c_ethanol_ag2,                                               &
        & hgrid, surface,                                                                                                &
        & t_cf('c_ethanol_ag2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_ag2', mem%c_nonsoluble_ag2,                                        &
        & hgrid, surface,                                                                                               &
        & t_cf('c_nonsoluble_ag2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_acid_bg2', mem%c_acid_bg2,                                               &
        & hgrid, surface,                                                                                          &
        & t_cf('c_acid_bg2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_bg2', mem%c_water_bg2,                                               &
        & hgrid, surface,                                                                                            &
        & t_cf('c_water_bg2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_bg2', mem%c_ethanol_bg2,                                              &
        & hgrid, surface,                                                                                               &
        & t_cf('c_ethanol_bg2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_bg2', mem%c_nonsoluble_bg2,                                       &
        & hgrid, surface,                                                                                              &
        & t_cf('c_nonsoluble_bg2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
        & prefix, suffix,                                                                                              &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_humus_2', mem%c_humus_2,                     &
        & hgrid, surface,                                                              &
        & t_cf('c_humus_2', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_bg_sum', mem%c_bg_sum,                                   &
        & hgrid, surface,                                                                &
        & t_cf('c_bg_sum', 'mol(C) m-2(canopy)', 'C-Pool of below ground C in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                                &
        & output_level=FULL,   loutput=.TRUE.,                                           &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL jsb_add_var( mem, 'c_ag_sum_1', mem%c_ag_sum_1,                                       &
        & hgrid, surface,                                                                        &
        & t_cf('c_ag_sum_1', 'mol(C) m-2(canopy)', 'C-Pool of above ground green C in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=FULL,  loutput=.TRUE.,                                                   &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

!      CALL jsb_add_var( mem, 'c_sum_natural', mem%C_sum_natural,                      &
!        & hgrid, surface,                                                             &
!        & t_cf('C_sum_natural', 'mol(C) m-2(canopy)', 'Sum of all natural C pools.'), &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
!        & prefix, suffix,                                                             &
!        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'lai_sum', mem%LAI_sum,                                   &
        & hgrid, surface,                                                              &
        & t_cf('LAI_sum', '-', 'Used to accumulate LAI over a day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_pot_sum', mem%NPP_pot_sum,                                   &
        & hgrid, surface,                                                                      &
        & t_cf('NPP_pot_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate NPP over a day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'gpp_sum', mem%GPP_sum,                                       &
        & hgrid, surface,                                                                  &
        & t_cf('GPP_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate GPP over a day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix,                                                                  &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

!@todo: relict from CBALANCE?
      CALL jsb_add_var( mem, 'topsoiltemp_sum', mem%topSoilTemp_sum,                                    &
        & hgrid, surface,                                                                               &
        & t_cf('topSoilTemp_sum', 'K', 'Used to accumulate upper soil layer temperature over a day..'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                            &
        & prefix, suffix,                                                                               &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

!@todo: relict from CBALANCE?
      CALL jsb_add_var( mem, 'alpha_sum', mem%alpha_sum,                                     &
        & hgrid, surface,                                                                    &
        & t_cf('alpha_sum', '-', 'Used to accumulate water stress coefficient over a day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'soil_respiration', mem%soil_respiration,           &
        & hgrid, surface,                                                        &
        & t_cf('soil_respiration', 'mol(C) m-2(canopy) s-1', 'Soil respiration.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_flux_correction', mem%NPP_flux_correction,                                      &
        & hgrid, surface,                                                                                         &
        & t_cf('NPP_flux_correction', 'mol(C) m-2(canopy) s-1', 'Amount by which the NPP rate has to be corrected &
        &                                                        to avoid a negative reserve pool.'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'excess_npp', mem%excess_NPP,                                                               &
        & hgrid, surface,                                                                                                &
        & t_cf('excess_NPP', 'mol(C) m-2 (canopy) s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'lai_yydaymean', mem%LAI_yyDayMean,              &
        & hgrid, surface,                                                     &
        & t_cf('LAI_yyDayMean', '-', 'Previous previous day mean.'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

     CALL jsb_add_var( mem, 'lai_ydaymean', mem%LAI_yDayMean,                 &
        & hgrid, surface,                                                     &
        & t_cf('LAI_yDayMean', '-', 'Previous day mean.'),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_pot_ydaymean', mem%NPP_pot_yDayMean,                              &
        & hgrid, surface,                                                                           &
        & t_cf('NPP_pot_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean NPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
        & prefix, suffix,                                                                           &
        & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                   ! in a zero divide in the first time step of an initialized run.

      CALL jsb_add_var( mem, 'npp_act_ydaymean', mem%NPP_act_yDayMean,                                &
        & hgrid, surface,                                                                             &
        & t_cf('NPP_act_yDayMean', 'mol(C) m-2(canopy) s-1', 'Actual NPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
        & prefix, suffix,                                                                             &
        & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                   ! in a zero divide in the first time step of an initialized run.

      CALL jsb_add_var( mem, 'gpp_ydaymean', mem%GPP_yDayMean,                                  &
        & hgrid, surface,                                                                       &
        & t_cf('GPP_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean GPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

!@todo: relict from CBALANCE?
      CALL jsb_add_var( mem, 'topsoiltemp_ydaymean', mem%topSoilTemp_yDayMean,                             &
        & hgrid, surface,                                                                                  &
        & t_cf('topSoilTemp_yDayMean', 'K', 'Previous day mean temperature of the uppermost soil layer.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                               &
        & prefix, suffix,                                                                                  &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

!@todo: relict from CBALANCE?
      CALL jsb_add_var( mem, 'alpha_ydaymean', mem%alpha_yDayMean,                             &
        & hgrid, surface,                                                                      &
        & t_cf('alpha_yDayMean', '-', 'Previous day mean value of Water Stress Coefficient.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix, loutput=.TRUE.,                                                                     &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3.

      CALL jsb_add_var( mem, 'root_exudates', mem%root_exudates,             &
        & hgrid, surface,                                                    &
        & t_cf('root_exudates', 'mol(C) m-2(canopy) s-1', 'Root exudates'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cflux_c_green_2_herb', mem%cflux_c_green_2_herb,                                              &
        & hgrid, surface,                                                                                         &
        & t_cf('cflux_c_green_2_herb', 'mol(C) m-2(canopy) s-1', 'Total carbon flux from C pool green to herbivory.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cflux_c_green_2_herb_lg', mem%cflux_herb_2_littergreen,                                         &
        & hgrid, surface,                                                                                          &
        & t_cf('cflux_herb_2_littergreen', 'mol(C) m-2(canopy) s-1', 'Carbon flux from cflux_c_green_2_herb into green litter &
        & (complementary to cflux_herb_2_atm).'),                                                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cflux_herb_2_atm', mem%cflux_herb_2_atm,                                    &
        & hgrid, surface,                                                                                           &
        & t_cf('cflux_herb_2_atm', 'mol(C) m-2(canopy) s-1', 'Carbon flux from cflux_c_green_2_herb into atmosphere &
        &  (complementary to cflux_herb_2_littergreen).'),                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

!       CALL jsb_add_var( mem, 'fract_litter_wood_new', mem%fract_litter_wood_new,                     &
!         & hgrid, surface,                                                                          &
!         & t_cf('fract_litter_wood_new', 'mol(C) m-2(canopy)', 'New fraction in wood litter pool.'), &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!         & prefix, suffix,                                                                          &
!         & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2greenlitterpools', mem%cflux_dist_greenreserve_2_soil,                           &
        & hgrid, surface,                                                                                          &
        & t_cf('cflux_dist_greenreserve_2_soil', 'mol(C) m-2(canopy) s-1', 'Amount of C relocated by wind damage to the &
        &  green litter pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                     &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2woodlitterpools', mem%cflux_dist_woods_2_soil,                             &
        & hgrid, surface,                                                                                          &
        & t_cf('cflux_dist_woods_2_soil', 'mol(C) m-2(canopy) s-1', 'Amount of C relocated by wind damage to the  &
        & woody litter pools..'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                     &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2atmos', mem%cflux_fire_all_2_atm,                                                  &
        & hgrid, surface,                                                                                           &
        & t_cf('cflux_fire_all_2_atm', 'mol(C) m-2(canopy) s-1', 'Carbon flux into atmosphere.'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=MEDIUM,                                                                                      &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3




      ! Carbon variables on tile area
      CALL jsb_add_var( mem, 'c_greenwood2litter_ta', mem%cflux_c_greenwood_2_litter_ta,                 &
        & hgrid, surface,                                                                                &
        & t_cf('c_greenwood2litter_ta', 'mol(CO2) m-2 s-1', 'Total litter flux entering the soil pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
        & prefix, suffix,                                                                                &
        & output_level=BASIC,                                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_green_ta', mem%c_green_ta,                              &
        & hgrid, surface,                                                                   &
        & t_cf('c_green_ta', 'mol(C) m-2', 'C-Pool for green parts of vegetation.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_woods_ta', mem%c_woods_ta,                                         &
        & hgrid, surface,                                                                              &
        & t_cf('c_woods_ta', 'mol(C) m-2', 'C-Pool for structural material/parts of vegetation.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput = .FALSE., output_level=BASIC,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_reserve_ta', mem%c_reserve_ta,                                                      &
        & hgrid, surface,                                                                                               &
        & t_cf('c_reserve_ta', 'mol(C) m-2', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=BASIC,                                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_crop_harvest_ta', mem%c_crop_harvest_ta,                        &
        & hgrid, surface,                                                                           &
        & t_cf('c_crop_harvest_ta', 'mol(C) m-2', 'C-Pool for Harvested Material from Crops.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
        & prefix, suffix,                                                                           &
        & output_level=BASIC,                                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_litter_grn_ag_ta', mem%c_litter_green_ag_ta,                              &
        & hgrid, surface,                                                                                       &
        & t_cf('c_litter_green_ag_ta', 'mol(C) m-2', 'C-Pool for above ground non-woody litter (leaves).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=BASIC, loutput=.FALSE.,                                                                 &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval in JSBACH3 is 60.0_wp

      CALL jsb_add_var( mem, 'c_litter_grn_bg_ta', mem%c_litter_green_bg_ta,                                  &
        & hgrid, surface,                                                                                           &
        & t_cf('c_litter_green_bg_ta', 'mol(C) m-2', 'C-Pool for below ground non-woody litter (fine roots).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC, loutput=.FALSE.,                                                                   &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval in JSBACH3 is 60.0_wp

      CALL jsb_add_var( mem, 'c_litter_wood_ag_ta', mem%c_litter_wood_ag_ta,                  &
        & hgrid, surface,                                                                         &
        & t_cf('c_litter_wood_ag_ta', 'mol(C) m-2', 'C-Pool for above ground woody litter.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      &
        & prefix, suffix,                                                                         &
        & output_level=BASIC, loutput=.FALSE.,                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval in JSBACH3 is 180.0_wp!

      CALL jsb_add_var( mem, 'c_litter_wood_bg_ta', mem%c_litter_wood_bg_ta,                  &
        & hgrid, surface,                                                                         &
        & t_cf('c_litter_wood_bg_ta', 'mol(C) m-2', 'C-Pool for below ground woody litter.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      &
        & prefix, suffix,                                                                         &
        & output_level=BASIC, loutput=.FALSE.,                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval in JSBACH3 is 180.0_wp

!CBAL       CALL jsb_add_var( mem, 'c_slow_ta', mem%c_slow_ta,                                          &
!CBAL         & hgrid, surface,                                                                             &
!CBAL         & t_cf('c_slow_ta', 'mol(C) m-2', 'C-Pool for slowly respirated soil organic material.'), &
!CBAL         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
!CBAL         & prefix, suffix,                                                                             &
!CBAL         & output_level=BASIC,                                                                            &
!CBAL         & lrestart=.FALSE., initval_r=2400.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3 from cbalance%c_slow.
!CBAL                                                  ! Note for cbalone this value would be:
!CBAL                                                  ! mo_cbalone_memory.f90:247: cbalance%c_slow(:,:) = 0.0_dp

      ! Yasso 1
      CALL jsb_add_var( mem, 'c_acid_ag1_ta', mem%c_acid_ag1_ta,                                         &
        & hgrid, surface,                                                                                     &
        & t_cf('c_acid_ag1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_ag1_ta', mem%c_water_ag1_ta,                                         &
        & hgrid, surface,                                                                                       &
        & t_cf('c_water_ag1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_ag1_ta', mem%c_ethanol_ag1_ta,                                         &
        & hgrid, surface,                                                                                           &
        & t_cf('c_ethanol_ag1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=MEDIUM,  loutput=.FALSE.,                                                                 &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_ag1_ta', mem%c_nonsoluble_ag1_ta,                                  &
        & hgrid, surface,                                                                                          &
        & t_cf('c_nonsoluble_ag1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_acid_bg1_ta', mem%c_acid_bg1_ta,                                         &
        & hgrid, surface,                                                                                     &
        & t_cf('c_acid_bg1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_bg1_ta', mem%c_water_bg1_ta,                                         &
        & hgrid, surface,                                                                                       &
        & t_cf('c_water_bg1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_bg1_ta', mem%c_ethanol_bg1_ta,                                        &
        & hgrid, surface,                                                                                          &
        & t_cf('c_ethanol_bg1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_bg1_ta', mem%c_nonsoluble_bg1_ta,                                 &
        & hgrid, surface,                                                                                         &
        & t_cf('c_nonsoluble_bg1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_humus_1_ta', mem%c_humus_1_ta,               &
        & hgrid, surface,                                                         &
        & t_cf('c_humus_1_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
        & prefix, suffix,                                                         &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      ! Yasso 2
      CALL jsb_add_var( mem, 'c_acid_ag2_ta', mem%c_acid_ag2_ta,                                         &
        & hgrid, surface,                                                                                     &
        & t_cf('c_acid_ag2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_ag2_ta', mem%c_water_ag2_ta,                                         &
        & hgrid, surface,                                                                                       &
        & t_cf('c_water_ag2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_ag2_ta', mem%c_ethanol_ag2_ta,                                         &
        & hgrid, surface,                                                                                           &
        & t_cf('c_ethanol_ag2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_ag2_ta', mem%c_nonsoluble_ag2_ta,                                  &
        & hgrid, surface,                                                                                          &
        & t_cf('c_nonsoluble_ag2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_acid_bg2_ta', mem%c_acid_bg2_ta,                                         &
        & hgrid, surface,                                                                                     &
        & t_cf('c_acid_bg2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_bg2_ta', mem%c_water_bg2_ta,                                         &
        & hgrid, surface,                                                                                       &
        & t_cf('c_water_bg2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_bg2_ta', mem%c_ethanol_bg2_ta,                                        &
        & hgrid, surface,                                                                                          &
        & t_cf('c_ethanol_bg2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_bg2_ta', mem%c_nonsoluble_bg2_ta,                                 &
        & hgrid, surface,                                                                                         &
        & t_cf('c_nonsoluble_bg2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_humus_2_ta', mem%c_humus_2_ta,               &
        & hgrid, surface,                                                         &
        & t_cf('c_humus_2_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
        & prefix, suffix,                                                         &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

!      CALL jsb_add_var( mem, 'c_bg_sum_ta', mem%c_bg_sum_ta,                         &
!        & hgrid, surface,                                                            &
!        & t_cf('c_bg_sum_ta', 'mol(C) m-2', 'C-Pool of below ground C in Yasso.'),   &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
!        & prefix, suffix,                                                            &
!        & output_level=MEDIUM,                                                       &
!        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL jsb_add_var( mem, 'c_sum_natural_ta', mem%C_sum_natural_ta,           &
        & hgrid, surface,                                                        &
        & t_cf('C_sum_natural_ta', 'mol(C) m-2', 'Sum of all natural C pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

!      CALL jsb_add_var( mem, 'npp_pot_sum_ta', mem%NPP_pot_sum_ta,                        &
!        & hgrid, surface,                                                                 &
!        & t_cf('NPP_pot_sum_ta', 'mol(C) m-2 s-1', 'Used to accumulate NPP over a day.'), &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
!        & prefix, suffix,                                                                 &
!        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

!      CALL jsb_add_var( mem, 'gpp_sum_ta', mem%GPP_sum_ta,                            &
!        & hgrid, surface,                                                             &
!        & t_cf('GPP_sum_ta', 'mol(C) m-2 s-1', 'Used to accumulate GPP over a day.'), &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
!        & prefix, suffix,                                                             &
!        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'soil_respiration_ta', mem%soil_respiration_ta,  &
        & hgrid, surface,                                                     &
        & t_cf('soil_respiration_ta', 'mol(C) m-2 s-1', 'Soil respiration.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_flux_corr_ta', mem%NPP_flux_correction_ta,                           &
        & hgrid, surface,                                                                                    &
        & t_cf('NPP_flux_correction_ta', 'mol(C) m-2 s-1', 'Amount by which the NPP rate has to be corrected &
        &                                                        to avoid a negative reserve pool.'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
        & prefix, suffix,                                                                                    &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'excess_npp_ta', mem%excess_NPP_ta,                                                    &
        & hgrid, surface,                                                                                           &
        & t_cf('excess_NPP_ta', 'mol(C) m-2  s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

!      CALL jsb_add_var( mem, 'lai_yydaymean_ta', mem%LAI_yyDayMean_ta,              &
!        & hgrid, surface,                                                     &
!        & t_cf('LAI_yyDayMean_ta', '-', 'Previous previous day mean.'),          &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
!        & prefix, suffix,                                                     &
!        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

!     CALL jsb_add_var( mem, 'lai_ydaymean_ta', mem%LAI_yDayMean_ta,                 &
!        & hgrid, surface,                                                     &
!        & t_cf('LAI_yDayMean_ta', '-', 'Previous day mean.'),                    &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
!        & prefix, suffix,                                                     &
!        & output_level=BASIC,                                                                                &
!        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_pot_ydaymean_ta', mem%NPP_pot_yDayMean_ta,                              &
        & hgrid, surface,                                                                           &
        & t_cf('NPP_pot_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean NPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
        & prefix, suffix,                                                                           &
        & output_level=BASIC, loutput = .TRUE.,                                                     &
        & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                   ! in a zero divide in the first time step of an initialized run.

      CALL jsb_add_var( mem, 'npp_act_ydaymean_ta', mem%NPP_act_yDayMean_ta,                                &
        & hgrid, surface,                                                                             &
        & t_cf('NPP_act_yDayMean_ta', 'mol(C) m-2 s-1', 'Actual NPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
        & prefix, suffix,                                                                             &
        & output_level=BASIC,loutput = .TRUE.,                                                          &
        & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                   ! in a zero divide in the first time step of an initialized run.

      CALL jsb_add_var( mem, 'gpp_ydaymean_ta', mem%GPP_yDayMean_ta,                       &
        & hgrid, surface,                                                                  &
        & t_cf('GPP_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean GPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix,                                                                  &
        & output_level=BASIC,loutput = .TRUE.,                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


      CALL jsb_add_var( mem, 'root_exudates_ta', mem%root_exudates_ta,       &
        & hgrid, surface,                                                    &
        & t_cf('root_exudates_ta', 'mol(C) m-2 s-1', 'Root exudates'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


      CALL jsb_add_var( mem, 'cflux_c_green_2_herb_ta', mem%cflux_c_green_2_herb_ta,                                   &
        & hgrid, surface,                                                                                    &
        & t_cf('cflux_c_green_2_herb_ta', 'mol(C) m-2 s-1', 'Total carbon flux from C pool green to herbivory.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
        & prefix, suffix,                                                                                    &
        & output_level=MEDIUM,                                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2greenlitter_ta', mem%cflux_dist_greenreserve_2_soil_ta,    &
        & hgrid, surface,                                                                         &
        & t_cf('cflux_dist_greenreserve_2_soil_ta', 'mol(C) m-2 s-1', 'Carbon flux to green litter.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      &
        & prefix, suffix,                                                                         &
        & output_level=MEDIUM,                                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2woodlitter_ta', mem%cflux_dist_woods_2_soil_ta,     &
        & hgrid, surface,                                                                        &
        & t_cf('cflux_dist_woods_2_soil_ta', 'mol(C) m-2 s-1', 'Carbon flux to woody litter.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=MEDIUM,                                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2atmos_ta', mem%cflux_fire_all_2_atm_ta,                                    &
        & hgrid, surface,                                                                                   &
        & t_cf('cflux_fire_all_2_atm_ta', 'mol(C) m-2 s-1', 'Total carbon flux from C pool green to herbivory.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
        & prefix, suffix,                                                                                   &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


!       CALL jsb_add_var( mem, 'cflux_herb_2_littergreen_ta', mem%cflux_herb_2_littergreen_ta,                                   &
!         & hgrid, surface,                                                                                          &
!         & t_cf('cflux_herb_2_littergreen_ta', 'mol(C) m-2 s-1', 'Carbon flux from cflux_c_green_2_herb into green litter &
!         &                                                      (complementary to cflux_herb_2_atm).'),        &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
!         & prefix, suffix,                                                                                          &
!         & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3
!
!       CALL jsb_add_var( mem, 'cflux_herb_2_atm_ta', mem%cflux_herb_2_atm_ta,                              &
!         & hgrid, surface,                                                                                           &
!         & t_cf('cflux_herb_2_atm_ta', 'mol(C) m-2 s-1', 'Carbon flux from cflux_c_green_2_herb into atmosphere &
!         &                                                          (complementary to cflux_herb_2_littergreen).'),        &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
!         & prefix, suffix,                                                                                           &
!         & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3


!       CALL jsb_add_var( mem, 'fract_litter_wood_new_ta', mem%fract_litter_wood_new_ta,                     &
!         & hgrid, surface,                                                                          &
!         & t_cf('fract_litter_wood_new_ta', 'mol(C) m-2', 'New fraction in wood litter pool.'), &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!         & prefix, suffix,                                                                          &
!         & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      ! R: In JSBACH3 the following variables are not written into the output.
      CALL jsb_add_var( mem, 'co2_l2a_npp_ta', mem%cflux_npp_2_atm_ta,                                                 &
        & hgrid, surface,                                                                                              &
        & t_cf('cflux_npp_2_atm_ta', 'kg(CO2) m-2 s-1', 'Net CO2 fluxes between biosphere and atmosphere due to NPP.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
        & prefix, suffix,                                                                                              &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'co2_l2a_resp_ta', mem%cflux_soilresp_2_atm_ta,                                 &
        & hgrid, surface,                                                                                    &
        & t_cf('cflux_soilresp_2_atm_ta', 'kg(CO2) m-2 s-1',                                                 &
        &      'Net CO2 fluxes between biosphere and atmosphere due to soil respiration.'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
        & prefix, suffix,                                                                                    &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'co2_l2a_herb_ta', mem%cflux_herb_2_atm_ta,                                 &
        & hgrid, surface,                                                                                     &
        & t_cf('cflux_herb_2_atm_ta', 'kg(CO2) m-2 s-1',                                                 &
        &      'Net CO2 fluxes between biosphere and atmosphere due to herbivory.'),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'co2_l2a_npp_yday_ta', mem%cflux_npp_2_atm_yday_ta,                                       &
        & hgrid, surface,                                                                                              &
        & t_cf('cflux_npp_2_atm_yday_ta', 'kg(CO2) m-2 s-1', 'Net daily CO2 flux, only needed for conservation test.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
        & prefix, suffix,                                                                                              &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ccons_calccpools', mem%cconservation_calcCpools,              &
        & hgrid, surface,                                                           &
        & t_cf('cconservation_calcCpools', 'mol(C) m-2(canopy)', &
          & 'Test for c conservation on allocation and respiration.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & loutput = .FALSE.,lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
          CALL jsb_add_var( mem,'max_green', mem%current_max_green,      &
            & hgrid, surface,                                                    &
            & t_cf('current_max_green', 'mol(C) m-2(canopy)', 'Running maximum of the green pool'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.)
          CALL jsb_add_var( mem,'max_green_veg_c', mem%veg_carbon_at_max_green, &
            & hgrid, surface,                                                    &
            & t_cf('veg_carbon_at_max_green', 'mol(C) m-2(canopy)', 'Total vegetation carbon at green pool maximum'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE.)
      ENDIF

    END IF


  END SUBROUTINE Init_carbon_memory

#endif
END MODULE mo_carbon_memory_class
