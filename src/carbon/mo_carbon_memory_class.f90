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
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE
  USE mo_jsb_process_class,      ONLY: FAGE_
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of prcesses in this module
  dsl4jsb_Use_processes FAGE_, PHENO_

  ! Use process configurations
  dsl4jsb_Use_config(PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_carbon_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 150

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
      litter_flux,                & ! R: From JS3: cbalance_diag_type. Carbon flux from the vegetation to the litter pools
                                    !    [mol(C)/m^2(canopy) s]
      Cpool_green,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts 
                                    ! .. of vegetation [mol(C)/m^2(canopy)]
      Cpool_woods,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                    ! .. plants [mol(C)/m^2(canopy)]
      Cpool_reserve,              & ! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive 
                                    ! .. bad times[mol(C)/m^2(canopy)]
!@todo why is the Cpool_crop_harvest realised as a density and not as a per tile variable in the first case?
      Cpool_crop_harvest,         & ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]     
     !SIZE CLASS 1
      YCpool_acid_ag1,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      YCpool_water_ag1,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      YCpool_ethanol_ag1,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      YCpool_nonsoluble_ag1,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      YCpool_acid_bg1,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      YCpool_water_bg1,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      YCpool_ethanol_bg1,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      YCpool_nonsoluble_bg1,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      YCpool_humus_1,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(canopy)]
     !SIZE CLASS 2
      YCpool_acid_ag2,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      YCpool_water_ag2,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      YCpool_ethanol_ag2,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      YCpool_nonsoluble_ag2,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      YCpool_acid_bg2,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
      YCpool_water_bg2,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
      YCpool_ethanol_bg2,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
      YCpool_nonsoluble_bg2,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
      YCpool_humus_2,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(canopy)]
      C_sum_natural,              & ! R: added
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
      Cflux_herbivory,            & ! Total carbon flux from C pool green to herbivory.
      Cflux_herbivory_LG,         & ! ..part of Cflux_herbivory that goes into green litter
      Cflux_herbivory_2_atm,      & ! ..part of Cflux_herbivory that goes into atmosphere
      !
      carbon_2_GreenLitterPools,  & ! R: should be later renamed to C_2_GreenLitterPools
      carbon_2_WoodLitterPools,   & ! R: should be later renamed to C_2_WoodLitterPools
      carbon_2_atmos,             &  ! R: should be later renamed to C_2_atmos
                                    ! Carbon immedieatly released by fire per tile [mol(C)m-2(subtile)]
      ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
      !frac_litter_wood_new,       & ! new fraction in above ground wood litter pool 
!JN: variables required when l_forestRegrowth = true, i.e. when the maxLAI is derived from biomass via allometric relationships
      current_max_green,          & ! Running maximum of the green pool [mol(C) m-2(canopy)] 
                                    ! since last maxLAI calculation
      veg_carbon_at_max_green       ! Total vegetation carbon [mol(C) m-2(canopy)] 
                                    ! at the time of the green pool maximum (current_max_green)

   ! Output on tile areas:
   TYPE(t_jsb_var_real2d) ::         &
!      litter_flux_ta,                & ! R: From JS3: cbalance_diag_type. Carbon flux from the vegetation to the litter pools
!                                       !    [mol(C)/m^2(tile area) s]
      Cpool_green_ta,                & ! C-pool for leaves, fine roots, vegetative organs and other green (living) parts 
                                       ! .. of vegetation [mol(C)/m^2(tile area)]
      Cpool_woods_ta,                & ! C-pool for stems, thick roots and other (dead) structural material of living
                                       ! .. plants [mol(C)/m^2(tile area)]
      Cpool_reserve_ta,              & ! C-pool for carbohydrate reserve (sugars, starches) that allows plants to survive 
                                       ! .. bad times[mol(C)/m^2(tile area)]
      Cpool_crop_harvest_ta,         & ! C-pool for biomass harvested from crops [mol(C)/m^2(grid box)]     
     !SIZE CLASS 1
      YCpool_acid_ag1_ta,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      YCpool_water_ag1_ta,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      YCpool_ethanol_ag1_ta,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      YCpool_nonsoluble_ag1_ta,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      YCpool_acid_bg1_ta,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      YCpool_water_bg1_ta,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      YCpool_ethanol_bg1_ta,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      YCpool_nonsoluble_bg1_ta,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      YCpool_humus_1_ta,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(tile area)]
     !SIZE CLASS 2
      YCpool_acid_ag2_ta,            & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      YCpool_water_ag2_ta,           & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      YCpool_ethanol_ag2_ta,         & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      YCpool_nonsoluble_ag2_ta,      & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      YCpool_acid_bg2_ta,            & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(tile area)]
      YCpool_water_bg2_ta,           & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(tile area)]
      YCpool_ethanol_bg2_ta,         & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(tile area)]
      YCpool_nonsoluble_bg2_ta,      & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(tile area)]
      YCpool_humus_2_ta,             & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(tile area)]
      C_sum_natural_ta,              & ! R: added
      !
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
      Cflux_herbivory_ta,            & ! Total carbon flux from C pool green to herbivory.
      !Cflux_herbivory_LG_ta,         & ! ..part of Cflux_herbivory that goes into green litter
      !Cflux_herbivory_2_atm_ta,      & ! ..part of Cflux_herbivory that goes into atmosphere
      carbon_2_GreenLitterPools_ta,  & ! R: should be later renamed to C_2_GreenLitterPools_ta
      carbon_2_WoodLitterPools_ta,   & ! R: should be later renamed to C_2_WoodLitterPools_ta
      carbon_2_atmos_ta,             & ! R: should be later renamed to C_2_atmos_ta
                                    ! Carbon immedieatly released by fire per tile [mol(C)m-2(subtile)]

      ! variables needed with Spitfire to compute fuel classes from wood and green pools and litter
      frac_litter_wood_new_ta,       & ! new fraction in above ground wood litter pool 
      ! former (JS3) local variables. R: nur auf tile area benötigt!
      CO2_flx2atm_npp_ta,            & ! [kg(CO2)/m^2 s]  tile averages of net CO2 fluxes between
      CO2_flx2atm_soilresp_ta,       & !  .. biosphere (due to NPP, soil respiration and
      CO2_flx2atm_herbivory_ta,      & !  .. grazing) and atmosphere [kg(CO2)/(m^2(ground) s)]
                                       ! .. (positive for emission to atmosphere,
                                       ! ..  negative for absorption by biosphere)
      CO2_flx2atm_npp_test_ta          ! CO2 flux due to NPP, only needed for conservation test. 

    ! Changes in state variables resulting due to area changes of the currently operating process 
    ! TODO-JN-COHORTS needs to be adapted for absolute fractions: currently in [mol(C)/m^2(parent)]
    ! TODO-JN currently only required for forests - but they do not have an lct_id
    TYPE(t_jsb_var_real2d) ::              &
      delta_Cpool_green_rp,                & 
      delta_Cpool_woods_rp,                & 
      delta_Cpool_reserve_rp,              & 
      delta_Cpool_crop_harvest_rp,         & 
      delta_YCpool_acid_ag1_rp,            & 
      delta_YCpool_water_ag1_rp,           & 
      delta_YCpool_ethanol_ag1_rp,         & 
      delta_YCpool_nonsoluble_ag1_rp,      & 
      delta_YCpool_acid_bg1_rp,            & 
      delta_YCpool_water_bg1_rp,           & 
      delta_YCpool_ethanol_bg1_rp,         & 
      delta_YCpool_nonsoluble_bg1_rp,      & 
      delta_YCpool_humus_1_rp,             & 
     !SIZE CLASS 2
      delta_YCpool_acid_ag2_rp,            & 
      delta_YCpool_water_ag2_rp,           & 
      delta_YCpool_ethanol_ag2_rp,         & 
      delta_YCpool_nonsoluble_ag2_rp,      & 
      delta_YCpool_acid_bg2_rp,            & 
      delta_YCpool_water_bg2_rp,           & 
      delta_YCpool_ethanol_bg2_rp,         & 
      delta_YCpool_nonsoluble_bg2_rp,      & 
      delta_YCpool_humus_2_rp,             &
      delta_current_max_green_rp,          & 
      delta_veg_carbon_at_max_green_rp,    &
      delta_carbon_2_GreenLitterPools_rp,  &
      delta_carbon_2_WoodLitterPools_rp,   &
      delta_carbon_2_atmos_rp


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
    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_carbon_memory'

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
        & lrestart=.FALSE., initval_r=0.0_wp )

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
        & lrestart=.TRUE., initval_r=3.0E-5_wp ) ! initval taken from JSBACH3, arbitary.

      ! R: In the output it should be the same as pseudo_precip, I do not understand why it was written out in JS3.
      CALL jsb_add_var( mem, 'pseudo_precip_yDay', mem%pseudo_precip_yDay,                             &
        & hgrid, surface,                                                                              &
        & t_cf('pseudo_precip_yDay', 'Kg m-2 s-1', 'Yesterdays pseudo-mean precipitation rate over previous 15 days'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & lrestart=.TRUE., initval_r=3.0E-5_wp ) ! initval taken from JSBACH3, arbitary.

      CALL jsb_add_var( mem, 'topsoiltemp', mem%topSoilTemp,                                           &
        & hgrid, surface,                                                                              &
        & t_cf('topSoilTemp', 'K', 'Top soil temperature.'),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=FULL,                                                                           &
        & lrestart=.TRUE., initval_r=283.15_wp ) ! I "invented" this initval.

      CALL jsb_add_var( mem, 'alpha', mem%alpha,                                                       &
        & hgrid, surface,                                                                              &
        & t_cf('alpha', '-', 'Water Stress Coefficient. Relative soil moisture.'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=FULL,                                                                           &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3's alpha_yDayMean


      ! Carbon variables on canopy area
      CALL jsb_add_var( mem, 'litter_flux', mem%litter_flux,                                             &
        & hgrid, surface,                                                                                &
        & t_cf('litter_flux', 'mol(CO2) m-2(canopy) s-1', 'Total litter flux entering the soil pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
        & prefix, suffix,                                                                                &
        & output_level=FULL,                                                                             &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cpool_green', mem%Cpool_green,                                     &
        & hgrid, surface,                                                                        &
        & t_cf('Cpool_green', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=FULL,                                                                             &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cpool_woods', mem%Cpool_woods,                                                &
        & hgrid, surface,                                                                                   &
        & t_cf('Cpool_woods', 'mol(C) m-2(canopy)', 'C-Pool for structural material/parts of vegetation.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
        & prefix, suffix,                                                                                   &
        & output_level=FULL,                                                                             &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cpool_reserve', mem%Cpool_reserve,                                                             &
        & hgrid, surface,                                                                                                    &
        & t_cf('Cpool_reserve', 'mol(C) m-2(canopy)', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                                 &
        & prefix, suffix,                                                                                                    &
        & output_level=FULL,                                                                             &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cpool_crop_harvest', mem%Cpool_crop_harvest,                               &
        & hgrid, surface,                                                                                &
        & t_cf('Cpool_crop_harvest', 'mol(C) m-2(canopy)', 'C-Pool for Harvested Material from Crops.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
        & prefix, suffix,                                                                                &
        & output_level=FULL,                                                                             &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      ! Yasso 1
      CALL jsb_add_var( mem, 'ycpool_acid_ag1', mem%YCpool_acid_ag1,                                               &
        & hgrid, surface,                                                                                          &
        & t_cf('YCpool_acid_ag1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=FULL,                                                                                       &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_water_ag1', mem%YCpool_water_ag1,                                               &
        & hgrid, surface,                                                                                            &
        & t_cf('YCpool_water_ag1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                                                       &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_ethanol_ag1', mem%YCpool_ethanol_ag1,                                               &
        & hgrid, surface,                                                                                                &
        & t_cf('YCpool_ethanol_ag1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & output_level=FULL,                                                                                       &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_nonsoluble_ag1', mem%YCpool_nonsoluble_ag1,                                        &
        & hgrid, surface,                                                                                               &
        & t_cf('YCpool_nonsoluble_ag1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=FULL,                                                                                       &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_acid_bg1', mem%YCpool_acid_bg1,                                               &
        & hgrid, surface,                                                                                          &
        & t_cf('YCpool_acid_bg1', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_water_bg1', mem%YCpool_water_bg1,                                               &
        & hgrid, surface,                                                                                            &
        & t_cf('YCpool_water_bg1', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_ethanol_bg1', mem%YCpool_ethanol_bg1,                                              &
        & hgrid, surface,                                                                                               &
        & t_cf('YCpool_ethanol_bg1', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_nonsoluble_bg1', mem%YCpool_nonsoluble_bg1,                                       &
        & hgrid, surface,                                                                                              &
        & t_cf('YCpool_nonsoluble_bg1', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
        & prefix, suffix,                                                                                              &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_humus_1', mem%YCpool_humus_1,                     &
        & hgrid, surface,                                                              &
        & t_cf('YCpool_humus_1', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      ! Yasso 2
      CALL jsb_add_var( mem, 'ycpool_acid_ag2', mem%YCpool_acid_ag2,                                               &
        & hgrid, surface,                                                                                          &
        & t_cf('YCpool_acid_ag2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_water_ag2', mem%YCpool_water_ag2,                                               &
        & hgrid, surface,                                                                                            &
        & t_cf('YCpool_water_ag2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_ethanol_ag2', mem%YCpool_ethanol_ag2,                                               &
        & hgrid, surface,                                                                                                &
        & t_cf('YCpool_ethanol_ag2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_nonsoluble_ag2', mem%YCpool_nonsoluble_ag2,                                        &
        & hgrid, surface,                                                                                               &
        & t_cf('YCpool_nonsoluble_ag2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_acid_bg2', mem%YCpool_acid_bg2,                                               &
        & hgrid, surface,                                                                                          &
        & t_cf('YCpool_acid_bg2', 'mol(C) m-2(canopy)', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_water_bg2', mem%YCpool_water_bg2,                                               &
        & hgrid, surface,                                                                                            &
        & t_cf('YCpool_water_bg2', 'mol(C) m-2(canopy)', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                         &
        & prefix, suffix,                                                                                            &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_ethanol_bg2', mem%YCpool_ethanol_bg2,                                              &
        & hgrid, surface,                                                                                               &
        & t_cf('YCpool_ethanol_bg2', 'mol(C) m-2(canopy)', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_nonsoluble_bg2', mem%YCpool_nonsoluble_bg2,                                       &
        & hgrid, surface,                                                                                              &
        & t_cf('YCpool_nonsoluble_bg2', 'mol(C) m-2(canopy)', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
        & prefix, suffix,                                                                                              &
        & output_level=FULL,                                                                                       &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ycpool_humus_2', mem%YCpool_humus_2,                     &
        & hgrid, surface,                                                              &
        & t_cf('YCpool_humus_2', 'mol(C) m-2(canopy)', 'C-Pool for humus in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & output_level=FULL,                                                                                       &
        & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_sum_natural', mem%C_sum_natural,                      &
        & hgrid, surface,                                                             &
        & t_cf('C_sum_natural', 'mol(C) m-2(canopy)', 'Sum of all natural C pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix,                                                             &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'lai_sum', mem%LAI_sum,                                   &
        & hgrid, surface,                                                              &
        & t_cf('LAI_sum', '-', 'Used to accumulate LAI over a day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
        & prefix, suffix,                                                              &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_pot_sum', mem%NPP_pot_sum,                                   &
        & hgrid, surface,                                                                      &
        & t_cf('NPP_pot_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate NPP over a day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'gpp_sum', mem%GPP_sum,                                       &
        & hgrid, surface,                                                                  &
        & t_cf('GPP_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate GPP over a day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix,                                                                  &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'topsoiltemp_sum', mem%topSoilTemp_sum,                                    &
        & hgrid, surface,                                                                               &
        & t_cf('topSoilTemp_sum', 'K', 'Used to accumulate upper soil layer temperature over a day..'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                            &
        & prefix, suffix,                                                                               &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

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
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_flux_correction', mem%NPP_flux_correction,                                      &
        & hgrid, surface,                                                                                         &
        & t_cf('NPP_flux_correction', 'mol(C) m-2(canopy) s-1', 'Amount by which the NPP rate has to be corrected &
        &                                                        to avoid a negative reserve pool.'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'excess_npp', mem%excess_NPP,                                                               &
        & hgrid, surface,                                                                                                &
        & t_cf('excess_NPP', 'mol(C) m-2 (canopy) s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

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
        & lrestart=.TRUE., initval_r=1.0E-13_wp ) ! initval taken from JSBACH3. Initializing to zero would result
                                                   ! in a zero divide in the first time step of an initialized run.

      CALL jsb_add_var( mem, 'npp_act_ydaymean', mem%NPP_act_yDayMean,                                &
        & hgrid, surface,                                                                             &
        & t_cf('NPP_act_yDayMean', 'mol(C) m-2(canopy) s-1', 'Actual NPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
        & prefix, suffix,                                                                             &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=1.0E-13_wp ) ! initval taken from JSBACH3. Initializing to zero would result
                                                   ! in a zero divide in the first time step of an initialized run.

      CALL jsb_add_var( mem, 'gpp_ydaymean', mem%GPP_yDayMean,                                  &
        & hgrid, surface,                                                                       &
        & t_cf('GPP_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean GPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
        & prefix, suffix,                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3


      CALL jsb_add_var( mem, 'topsoiltemp_ydaymean', mem%topSoilTemp_yDayMean,                             &
        & hgrid, surface,                                                                                  &
        & t_cf('topSoilTemp_yDayMean', 'K', 'Previous day mean temperature of the uppermost soil layer.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                               &
        & prefix, suffix,                                                                                  &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'alpha_ydaymean', mem%alpha_yDayMean,                             &
        & hgrid, surface,                                                                      &
        & t_cf('alpha_yDayMean', '-', 'Previous day mean value of Water Stress Coefficient.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                   &
        & prefix, suffix,                                                                      &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3.

      CALL jsb_add_var( mem, 'root_exudates', mem%root_exudates,             &
        & hgrid, surface,                                                    &
        & t_cf('root_exudates', 'mol(C) m-2(canopy) s-1', 'Root exudates'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cflux_herbivory', mem%Cflux_herbivory,                                              &
        & hgrid, surface,                                                                                         &
        & t_cf('Cflux_herbivory', 'mol(C) m-2(canopy) s-1', 'Total carbon flux from C pool green to herbivory.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cflux_herbivory_lg', mem%Cflux_herbivory_LG,                                         &
        & hgrid, surface,                                                                                          &
        & t_cf('Cflux_herbivory_LG', 'mol(C) m-2(canopy) s-1', 'Carbon flux from Cflux_herbivory into green litter &
        & (complementary to Cflux_herbivory_2_atm).'),                                                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'cflux_herbivory_2_atm', mem%Cflux_herbivory_2_atm,                                    &
        & hgrid, surface,                                                                                           &
        & t_cf('Cflux_herbivory_2_atm', 'mol(C) m-2(canopy) s-1', 'Carbon flux from Cflux_herbivory into atmosphere &
        &  (complementary to Cflux_herbivory_LG).'),                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

!       CALL jsb_add_var( mem, 'frac_litter_wood_new', mem%frac_litter_wood_new,                     &
!         & hgrid, surface,                                                                          &
!         & t_cf('frac_litter_wood_new', 'mol(C) m-2(canopy)', 'New fraction in wood litter pool.'), &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!         & prefix, suffix,                                                                          &
!         & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2greenlitterpools', mem%carbon_2_GreenLitterPools,                           &
        & hgrid, surface,                                                                                          &
        & t_cf('carbon_2_GreenLitterPools', 'mol(C) m-2(canopy) s-1', 'Amount of C relocated by wind damage to the &
        &  green litter pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                     &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2woodlitterpools', mem%carbon_2_WoodLitterPools,                             &
        & hgrid, surface,                                                                                          &
        & t_cf('carbon_2_WoodLitterPools', 'mol(C) m-2(canopy) s-1', 'Amount of C relocated by wind damage to the  &
        & woody litter pools..'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                     &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2atmos', mem%carbon_2_atmos,                                                  &
        & hgrid, surface,                                                                                           &
        & t_cf('carbon_2_atmos', 'mol(C) m-2(canopy) s-1', 'Carbon flux into atmosphere.'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=MEDIUM,                                                                                      &
        & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3



      ! Carbon variables on tile area
      ! CALL jsb_add_var( mem, 'litter_flux_ta', mem%litter_flux_ta,                                       &
      !   & hgrid, surface,                                                                                &
      !   & t_cf('litter_flux_ta', 'mol(CO2) m-2 s-1', 'Total litter flux entering the soil pools.'),      &
      !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      !   & prefix, suffix,                                                                                &
      !   & output_level=BASIC,                                                                            &
      !   & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_green_ta', mem%Cpool_green_ta,                              &
        & hgrid, surface,                                                                   &
        & t_cf('Cpool_green_ta', 'mol(C) m-2', 'C-Pool for green parts of vegetation.'),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                                            &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_woods_ta', mem%Cpool_woods_ta,                                         &
        & hgrid, surface,                                                                              &
        & t_cf('Cpool_woods_ta', 'mol(C) m-2', 'C-Pool for structural material/parts of vegetation.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level=BASIC,                                                                            &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_reserve_ta', mem%Cpool_reserve_ta,                                                      &
        & hgrid, surface,                                                                                               &
        & t_cf('Cpool_reserve_ta', 'mol(C) m-2', 'C-Pool for reserve carbohydrates (starches, sugars) of vegetation.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                            &
        & prefix, suffix,                                                                                               &
        & output_level=BASIC,                                                                            &
        & loutput = .TRUE.,lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_crop_harvest_ta', mem%Cpool_crop_harvest_ta,                        &
        & hgrid, surface,                                                                           &
        & t_cf('Cpool_crop_harvest_ta', 'mol(C) m-2', 'C-Pool for Harvested Material from Crops.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
        & prefix, suffix,                                                                           &
        & output_level=BASIC,                                                                            &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      ! Yasso 1
      CALL jsb_add_var( mem, 'c_acid_ag1_ta', mem%YCpool_acid_ag1_ta,                                         &
        & hgrid, surface,                                                                                     &
        & t_cf('YCpool_acid_ag1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=MEDIUM,                                                                                &
        & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_ag1_ta', mem%YCpool_water_ag1_ta,                                         &
        & hgrid, surface,                                                                                       &
        & t_cf('YCpool_water_ag1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_ag1_ta', mem%YCpool_ethanol_ag1_ta,                                         &
        & hgrid, surface,                                                                                           &
        & t_cf('YCpool_ethanol_ag1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_ag1_ta', mem%YCpool_nonsoluble_ag1_ta,                                  &
        & hgrid, surface,                                                                                          &
        & t_cf('YCpool_nonsoluble_ag1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_acid_bg1_ta', mem%YCpool_acid_bg1_ta,                                         &
        & hgrid, surface,                                                                                     &
        & t_cf('YCpool_acid_bg1_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_bg1_ta', mem%YCpool_water_bg1_ta,                                         &
        & hgrid, surface,                                                                                       &
        & t_cf('YCpool_water_bg1_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_bg1_ta', mem%YCpool_ethanol_bg1_ta,                                        &
        & hgrid, surface,                                                                                          &
        & t_cf('YCpool_ethanol_bg1_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_bg1_ta', mem%YCpool_nonsoluble_bg1_ta,                                 &
        & hgrid, surface,                                                                                         &
        & t_cf('YCpool_nonsoluble_bg1_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_humus_1_ta', mem%YCpool_humus_1_ta,               &
        & hgrid, surface,                                                         &
        & t_cf('YCpool_humus_1_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
        & prefix, suffix,                                                         &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      ! Yasso 2
      CALL jsb_add_var( mem, 'c_acid_ag2_ta', mem%YCpool_acid_ag2_ta,                                         &
        & hgrid, surface,                                                                                     &
        & t_cf('YCpool_acid_ag2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_ag2_ta', mem%YCpool_water_ag2_ta,                                         &
        & hgrid, surface,                                                                                       &
        & t_cf('YCpool_water_ag2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_ag2_ta', mem%YCpool_ethanol_ag2_ta,                                         &
        & hgrid, surface,                                                                                           &
        & t_cf('YCpool_ethanol_ag2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (aboveground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_ag2_ta', mem%YCpool_nonsoluble_ag2_ta,                                  &
        & hgrid, surface,                                                                                          &
        & t_cf('YCpool_nonsoluble_ag2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (aboveground).'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_acid_bg2_ta', mem%YCpool_acid_bg2_ta,                                         &
        & hgrid, surface,                                                                                     &
        & t_cf('YCpool_acid_bg2_ta', 'mol(C) m-2', 'C-Pool for acid-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_water_bg2_ta', mem%YCpool_water_bg2_ta,                                         &
        & hgrid, surface,                                                                                       &
        & t_cf('YCpool_water_bg2_ta', 'mol(C) m-2', 'C-Pool for water-soluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                    &
        & prefix, suffix,                                                                                       &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_ethanol_bg2_ta', mem%YCpool_ethanol_bg2_ta,                                        &
        & hgrid, surface,                                                                                          &
        & t_cf('YCpool_ethanol_bg2_ta', 'mol(C) m-2', 'C-Pool for ethanol-soluble litter in Yasso (belowground)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
        & prefix, suffix,                                                                                          &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_nonsoluble_bg2_ta', mem%YCpool_nonsoluble_bg2_ta,                                 &
        & hgrid, surface,                                                                                         &
        & t_cf('YCpool_nonsoluble_bg2_ta', 'mol(C) m-2', 'C-Pool for nonsoluble litter in Yasso (belowground).'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                      &
        & prefix, suffix,                                                                                         &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_humus_2_ta', mem%YCpool_humus_2_ta,               &
        & hgrid, surface,                                                         &
        & t_cf('YCpool_humus_2_ta', 'mol(C) m-2', 'C-Pool for humus in Yasso.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
        & prefix, suffix,                                                         &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c_sum_natural_ta', mem%C_sum_natural_ta,           &
        & hgrid, surface,                                                        &
        & t_cf('C_sum_natural_ta', 'mol(C) m-2', 'Sum of all natural C pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & output_level=MEDIUM,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      ! CALL jsb_add_var( mem, 'npp_pot_sum_ta', mem%NPP_pot_sum_ta,                        &
      !   & hgrid, surface,                                                                 &
      !   & t_cf('NPP_pot_sum_ta', 'mol(C) m-2 s-1', 'Used to accumulate NPP over a day.'), &
      !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
      !   & prefix, suffix,                                                                 &
      !   & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      ! CALL jsb_add_var( mem, 'gpp_sum_ta', mem%GPP_sum_ta,                            &
      !   & hgrid, surface,                                                             &
      !   & t_cf('GPP_sum_ta', 'mol(C) m-2 s-1', 'Used to accumulate GPP over a day.'), &
      !   & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
      !   & prefix, suffix,                                                             &
      !   & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'soil_respiration_ta', mem%soil_respiration_ta,  &
        & hgrid, surface,                                                     &
        & t_cf('soil_respiration_ta', 'mol(C) m-2 s-1', 'Soil respiration.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & output_level=BASIC,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_flux_corr_ta', mem%NPP_flux_correction_ta,                           &
        & hgrid, surface,                                                                                    &
        & t_cf('NPP_flux_correction_ta', 'mol(C) m-2 s-1', 'Amount by which the NPP rate has to be corrected &
        &                                                        to avoid a negative reserve pool.'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
        & prefix, suffix,                                                                                    &
        & output_level=BASIC,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'excess_npp_ta', mem%excess_NPP_ta,                                                    &
        & hgrid, surface,                                                                                           &
        & t_cf('excess_NPP_ta', 'mol(C) m-2  s-1', 'NPP that could not be not stored in vegetation carbon pools.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

     !  CALL jsb_add_var( mem, 'lai_yydaymean_ta', mem%LAI_yyDayMean_ta,              &
     !    & hgrid, surface,                                                     &
     !    & t_cf('LAI_yyDayMean_ta', '-', 'Previous previous day mean.'),          &
     !    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
     !    & prefix, suffix,                                                     &
     !    & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

     ! CALL jsb_add_var( mem, 'lai_ydaymean_ta', mem%LAI_yDayMean_ta,                 &
     !    & hgrid, surface,                                                     &
     !    & t_cf('LAI_yDayMean_ta', '-', 'Previous day mean.'),                    &
     !    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
     !    & prefix, suffix,                                                     &
     !    & output_level=BASIC,                                                                                &
     !    & lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_pot_ydaymean_ta', mem%NPP_pot_yDayMean_ta,                              &
        & hgrid, surface,                                                                           &
        & t_cf('NPP_pot_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean NPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                        &
        & prefix, suffix,                                                                           &
        & output_level=BASIC,                                                                                &
        & lrestart=.TRUE., initval_r=1.0E-13_wp ) ! initval taken from JSBACH3. Initializing to zero would result
                                                   ! in a zero divide in the first time step of an initialized run.

      CALL jsb_add_var( mem, 'npp_act_ydaymean_ta', mem%NPP_act_yDayMean_ta,                                &
        & hgrid, surface,                                                                             &
        & t_cf('NPP_act_yDayMean_ta', 'mol(C) m-2 s-1', 'Actual NPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
        & prefix, suffix,                                                                             &
        & output_level=BASIC,                                                                                &
        & lrestart=.TRUE., initval_r=1.0E-13_wp ) ! initval taken from JSBACH3. Initializing to zero would result
                                                   ! in a zero divide in the first time step of an initialized run.

      CALL jsb_add_var( mem, 'gpp_ydaymean_ta', mem%GPP_yDayMean_ta,                       &
        & hgrid, surface,                                                                  &
        & t_cf('GPP_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean GPP Rate of the previous day.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix,                                                                  &
        & loutput = .TRUE., &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3


      CALL jsb_add_var( mem, 'root_exudates_ta', mem%root_exudates_ta,       &
        & hgrid, surface,                                                    &
        & t_cf('root_exudates_ta', 'mol(C) m-2 s-1', 'Root exudates'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3


      CALL jsb_add_var( mem, 'cflux_herbivory_ta', mem%Cflux_herbivory_ta,                                   &
        & hgrid, surface,                                                                                    &
        & t_cf('Cflux_herbivory_ta', 'mol(C) m-2 s-1', 'Total carbon flux from C pool green to herbivory.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
        & prefix, suffix,                                                                                    &
        & output_level=MEDIUM,                                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2greenlitter_ta', mem%carbon_2_GreenLitterPools_ta,    &
        & hgrid, surface,                                                                         &
        & t_cf('carbon_2_GreenLitterPools_ta', 'mol(C) m-2 s-1', 'Carbon flux to green litter.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      &
        & prefix, suffix,                                                                         &
        & output_level=MEDIUM,                                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2woodlitter_ta', mem%carbon_2_WoodLitterPools_ta,     &
        & hgrid, surface,                                                                        &
        & t_cf('carbon_2_WoodLitterPools_ta', 'mol(C) m-2 s-1', 'Carbon flux to woody litter.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                                        &
        & output_level=MEDIUM,                                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'c2atmos_ta', mem%carbon_2_atmos_ta,                                    &
        & hgrid, surface,                                                                                   &
        & t_cf('carbon_2_atmos_ta', 'mol(C) m-2 s-1', 'Total carbon flux from C pool green to herbivory.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                &
        & prefix, suffix,                                                                                   &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3


!       CALL jsb_add_var( mem, 'Cflux_herbivory_LG_ta', mem%Cflux_herbivory_LG_ta,                                   &
!         & hgrid, surface,                                                                                          &
!         & t_cf('Cflux_herbivory_LG_ta', 'mol(C) m-2 s-1', 'Carbon flux from Cflux_herbivory into green litter &
!         &                                                      (complementary to Cflux_herbivory_2_atm).'),        &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                       &
!         & prefix, suffix,                                                                                          &
!         & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3
! 
!       CALL jsb_add_var( mem, 'Cflux_herbivory_2_atm_ta', mem%Cflux_herbivory_2_atm_ta,                              &
!         & hgrid, surface,                                                                                           &
!         & t_cf('Cflux_herbivory_2_atm_ta', 'mol(C) m-2 s-1', 'Carbon flux from Cflux_herbivory into atmosphere &
!         &                                                          (complementary to Cflux_herbivory_LG).'),        &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
!         & prefix, suffix,                                                                                           &
!         & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3


!       CALL jsb_add_var( mem, 'frac_litter_wood_new_ta', mem%frac_litter_wood_new_ta,                     &
!         & hgrid, surface,                                                                          &
!         & t_cf('frac_litter_wood_new_ta', 'mol(C) m-2', 'New fraction in wood litter pool.'), &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
!         & prefix, suffix,                                                                          &
!         & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      ! R: In JSBACH3 the following variables are not written into the output.
      CALL jsb_add_var( mem, 'co2_l2a_npp_ta', mem%CO2_flx2atm_npp_ta,                                                 &
        & hgrid, surface,                                                                                              &
        & t_cf('CO2_flx2atm_npp_ta', 'kg(CO2) m-2 s-1', 'Net CO2 fluxes between biosphere and atmosphere due to NPP.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                           &
        & prefix, suffix,                                                                                              &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'co2_l2a_resp_ta', mem%CO2_flx2atm_soilresp_ta,                                 &
        & hgrid, surface,                                                                                    &
        & t_cf('CO2_flx2atm_soilresp_ta', 'kg(CO2) m-2 s-1',                                                 &
        &      'Net CO2 fluxes between biosphere and atmosphere due to soil respiration.'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                 &
        & prefix, suffix,                                                                                    &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'co2_l2a_herb_ta', mem%CO2_flx2atm_herbivory_ta,                                 &
        & hgrid, surface,                                                                                     &
        & t_cf('CO2_flx2atm_herbivory_ta', 'kg(CO2) m-2 s-1',                                                 &
        &      'Net CO2 fluxes between biosphere and atmosphere due to herbivory.'),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                  &
        & prefix, suffix,                                                                                     &
        & output_level=BASIC,                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

!       CALL jsb_add_var( mem, 'CO2_flx2atm_npp_test_ta', mem%CO2_flx2atm_npp_test_ta,                                           &
!         & hgrid, surface,                                                                                                      &
!         & t_cf('CO2_flx2atm_npp_test_ta', 'mol(C) m-2 s-1', 'CO2 flux due to NPP, only needed for conservation test'), &
!         & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                                   &
!         & prefix, suffix,                                                                                                      &
!         & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'ccons_calccpools', mem%cconservation_calcCpools,              &
        & hgrid, surface,                                                           &
        & t_cf('cconservation_calcCpools', 'mol(C) m-2(canopy)', &
          & 'Test for c conservation on allocation and respiration.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp ) 

      IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
          CALL jsb_add_var( mem,'max_green', mem%current_max_green,      &
            & hgrid, surface,                                                    &
            & t_cf('current_max_green', 'mol(C) m-2(canopy)', 'Running maximum of the green pool'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)
          CALL jsb_add_var( mem,'max_green_veg_c', mem%veg_carbon_at_max_green, &
            & hgrid, surface,                                                    &
            & t_cf('veg_carbon_at_max_green', 'mol(C) m-2(canopy)', 'Total vegetation carbon at green pool maximum'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp)
      ENDIF

      IF (model%processes(FAGE_)%p%config%active) THEN

        IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
          CALL jsb_add_var( mem, 'delta_current_max_green_rp', mem%delta_current_max_green_rp, hgrid, surface, &
            & t_cf('delta_current_max_green_rp', 'mol(C) m-2(parent))', 'Change in current max green.'),         &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
            & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )

          CALL jsb_add_var( mem, 'delta_veg_carbon_at_max_green_rp', mem%delta_veg_carbon_at_max_green_rp, hgrid, surface, &
            & t_cf('delta_veg_carbon_at_max_green_rp', 'mol(C) m-2(parent))', 'Change in veg at current max green.'),      &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
            & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        ENDIF


!TODO-JN-COHORTS: does not need to go to the output (loutput = .FALSE.) -- only now for testing
!TODO-JN-COHORTS: change unit when absolute areas are implemented
! ++ think about "fract_fpc_max"!
        CALL jsb_add_var( mem, 'delta_cpool_green_rp', mem%delta_Cpool_green_rp, hgrid, surface, &
          & t_cf('delta_Cpool_green_rp', 'mol(C) m-2(parent))', 'Change in Cpool_green.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_cpool_woods_rp', mem%delta_Cpool_woods_rp, hgrid, surface, &
          & t_cf('delta_Cpool_woods_rp', 'mol(C) m-2(parent))', 'Change in Cpool_woods.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_cpool_reserve_rp', mem%delta_Cpool_reserve_rp, hgrid, surface, &
          & t_cf('delta_Cpool_reserve_rp', 'mol(C) m-2(parent))', 'Change in Cpool_reserve.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp )

        CALL jsb_add_var( mem, 'delta_cpool_crop_harvest_rp', mem%delta_Cpool_crop_harvest_rp, hgrid, surface, &
          & t_cf('delta_Cpool_crop_harvest_rp', 'mol(C) m-2(parent))', 'Change in Cpool_crop_harvest.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )

        CALL jsb_add_var( mem, 'delta_ycpool_acid_ag1_rp', mem%delta_YCpool_acid_ag1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_acid_ag1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_acid_ag1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_water_ag1_rp', mem%delta_YCpool_water_ag1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_water_ag1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_water_ag1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_ethanol_ag1_rp', mem%delta_YCpool_ethanol_ag1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_ethanol_ag1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_ethanol_ag1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_nonsoluble_ag1_rp', mem%delta_YCpool_nonsoluble_ag1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_nonsoluble_ag1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_nonsoluble_ag1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_acid_ag2_rp', mem%delta_YCpool_acid_ag2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_acid_ag2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_acid_ag2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_water_ag2_rp', mem%delta_YCpool_water_ag2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_water_ag2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_water_ag2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_ethanol_ag2_rp', mem%delta_YCpool_ethanol_ag2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_ethanol_ag2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_ethanol_ag2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_nonsoluble_ag2_rp', mem%delta_YCpool_nonsoluble_ag2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_nonsoluble_ag2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_nonsoluble_ag2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_humus_1_rp', mem%delta_YCpool_humus_1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_humus_1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_humus_1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )

        CALL jsb_add_var( mem, 'delta_ycpool_acid_bg1_rp', mem%delta_YCpool_acid_bg1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_acid_bg1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_acid_bg1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_water_bg1_rp', mem%delta_YCpool_water_bg1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_water_bg1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_water_bg1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_ethanol_bg1_rp', mem%delta_YCpool_ethanol_bg1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_ethanol_bg1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_ethanol_bg1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_YCpool_nonsoluble_bg1_rp', mem%delta_YCpool_nonsoluble_bg1_rp, hgrid, surface, &
          & t_cf('delta_YCpool_nonsoluble_bg1_rp', 'mol(C) m-2(parent))', 'Change in YCpool_nonsoluble_bg1.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_acid_bg2_rp', mem%delta_YCpool_acid_bg2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_acid_bg2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_acid_bg2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_water_bg2_rp', mem%delta_YCpool_water_bg2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_water_bg2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_water_bg2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_ethanol_bg2_rp', mem%delta_YCpool_ethanol_bg2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_ethanol_bg2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_ethanol_bg2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_nonsoluble_bg2_rp', mem%delta_YCpool_nonsoluble_bg2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_nonsoluble_bg2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_nonsoluble_bg2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_ycpool_humus_2_rp', mem%delta_YCpool_humus_2_rp, hgrid, surface, &
          & t_cf('delta_YCpool_humus_2_rp', 'mol(C) m-2(parent))', 'Change in YCpool_humus_2.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )

        CALL jsb_add_var( mem, 'delta_c2glpools_rp', mem%delta_carbon_2_GreenLitterPools_rp, hgrid, surface, &
          & t_cf('delta_carbon_2_GreenLitterPools_rp', 'mol(C) m-2(parent))', 'Change in carbon_2_GreenLitterPools.'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_c2wlpools_rp', mem%delta_carbon_2_WoodLitterPools_rp, hgrid, surface, &
          & t_cf('delta_carbon_2_WoodLitterPools_rp', 'mol(C) m-2(parent))', 'Change in carbon_2_WoodLitterPools.'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )
        CALL jsb_add_var( mem, 'delta_c2atmos_rp', mem%delta_carbon_2_atmos_rp, hgrid, surface, &
          & t_cf('delta_carbon_2_atmos_rp', 'mol(C) m-2(parent))', 'Change in carbon_2_atmos.'),         &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
          & prefix, suffix, loutput = .FALSE., lrestart=.FALSE., initval_r=0.0_wp )


      END IF ! (model%processes(FAGE_)%p%config%active) THEN

    END IF


  END SUBROUTINE Init_carbon_memory

#endif
END MODULE mo_carbon_memory_class
