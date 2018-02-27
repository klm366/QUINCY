!>
!! @brief Contains the routines for the carbon processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!! Julia Nabel (2018-01-10): l_forestRegrowth flag -> allow for changing max_LAI
!!
MODULE mo_carbon_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_lctlib_class, ONLY: t_lctlib_element

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_Cpools, relocate_carbon_fire, relocate_carbon_damage, yasso, &
    & update_delta_carbon, update_state_reset_delta, get_per_tile

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_process'

CONTAINS

  ! Calculations to allocate NPP to the plant C pools.
  elemental subroutine calc_Cpools( LAI,                       & ! in
                                    LAI_previous,              & ! in
                                    max_LAI,                   & ! in
                                    NPP_pot_yDayMean,          & ! in
                                    topSoilTemp,               & ! in
                                    alpha_yDayMean,            & ! in
                                    !
                                    frac_npp_2_woodPool,       & ! in
                                    frac_npp_2_reservePool,    & ! in
                                    frac_npp_2_exudates,       & ! in
                                    frac_green_2_herbivory,    & ! in
                                    tau_Cpool_litter_green,    & ! in
                                    tau_Cpool_litter_wood,     & ! in
                                    tau_Cpool_woods,           & ! in
                                    LAI_shed_constant,         & ! in
                                    frac_C_litter_green2atmos, & ! in
                                    Max_C_content_woods,       & ! in
                                    specific_leaf_area_C,      & ! in
                                    reserveC2leafC,            & ! in
                                    !is_vegetation,             & ! in
                                    pheno_type,                & ! in
                                    time_steps_per_day,        & ! in
                                    dtime,                     & ! in
                                    !
                                    Cpool_green,               & ! inout
                                    Cpool_woods,               & ! inout
                                    Cpool_reserve,             & ! inout
                                    Cpool_crop_harvest,        & ! inout
                                    !
                                    soilResp_rate,             & ! out
                                    NPP_flux_correction,       & ! out
                                    excess_NPP,                & ! out
                                    root_exudates,             & ! out
                                    Cflx_litterTotal,          & ! out
                                    Cflx_herbivory,            & ! out
                                    Cflx_herbivory_LG,         & ! out
                                    Cflx_herbivory_2_atm,      & ! out
                                    NPP_act_yDayMean,          & ! out
                                    !
                                    !frac_litter_wood_new,      & ! inout, optional; R: only for spitfire
                                    ! Yasso variables
                                    temp2_30d,                 & ! in, optional, necessary with yasso
                                    precip_30d,                & ! in, optional, necessary with yasso
                                    !
                                    YCpool_acid_ag1,           & ! inout, optional, necessary with yasso
                                    YCpool_water_ag1,          & ! inout, optional, necessary with yasso
                                    YCpool_ethanol_ag1,        & ! inout, optional, necessary with yasso
                                    YCpool_nonsoluble_ag1,     & ! inout, optional, necessary with yasso
                                    YCpool_acid_bg1,           & ! inout, optional, necessary with yasso
                                    YCpool_water_bg1,          & ! inout, optional, necessary with yasso
                                    YCpool_ethanol_bg1,        & ! inout, optional, necessary with yasso
                                    YCpool_nonsoluble_bg1,     & ! inout, optional, necessary with yasso
                                    YCpool_humus_1,            & ! inout, optional, necessary with yasso
                                    YCpool_acid_ag2,           & ! inout, optional, necessary with yasso
                                    YCpool_water_ag2,          & ! inout, optional, necessary with yasso
                                    YCpool_ethanol_ag2,        & ! inout, optional, necessary with yasso
                                    YCpool_nonsoluble_ag2,     & ! inout, optional, necessary with yasso
                                    YCpool_acid_bg2,           & ! inout, optional, necessary with yasso
                                    YCpool_water_bg2,          & ! inout, optional, necessary with yasso
                                    YCpool_ethanol_bg2,        & ! inout, optional, necessary with yasso
                                    YCpool_nonsoluble_bg2,     & ! inout, optional, necessary with yasso
                                    YCpool_humus_2,            & ! inout, optional, necessary with yasso
                                    LeafLit_coef1,             & ! in, optional, necessary with yasso
                                    LeafLit_coef2,             & ! in, optional, necessary with yasso
                                    LeafLit_coef3,             & ! in, optional, necessary with yasso
                                    LeafLit_coef4,             & ! in, optional, necessary with yasso
                                    LeafLit_coef5,             & ! in, optional, necessary with yasso
                                    WoodLit_coef1,             & ! in, optional, necessary with yasso
                                    WoodLit_coef2,             & ! in, optional, necessary with yasso
                                    WoodLit_coef3,             & ! in, optional, necessary with yasso
                                    WoodLit_coef4,             & ! in, optional, necessary with yasso
                                    WoodLit_coef5,             & ! in, optional, necessary with yasso
                                    WoodLitterSize             & ! in, optional, necessary with yasso
                                    )

    USE mo_carbon_constants,       ONLY : frac_green_aboveGround, frac_wood_aboveGround
    USE mo_jsb_physical_constants, ONLY : tmelt     ! tmelt =273.15_wp
    USE mo_carbon_constants,       ONLY : sec_per_day, days_per_year, sec_per_year

    ! Arguments
    real(wp),intent(in)    :: LAI                    ! Yesterdays mean LAI
    real(wp),intent(in)    :: LAI_previous           ! The day before yesterdays mean LAI
    real(wp),intent(in)    :: max_LAI                ! Maximum value of LAI
    real(wp),intent(in)    :: NPP_pot_yDayMean       ! Yesterdays mean NPP rate [mol(C)/m^2 s]
    real(wp),intent(in)    :: topSoilTemp            ! Yesterdays mean temperature of upper soil layer [degree Kelvin]
    real(wp),intent(in)    :: alpha_yDayMean         ! Yesterdays mean water stress factor (between 0 and 1)
    real(wp),intent(in)    :: frac_npp_2_woodPool    ! Fraction of NPP to be put maximally into the green pool
    real(wp),intent(in)    :: frac_npp_2_reservePool ! Fraction of NPP to be put into the optimally into the reserve pool
    real(wp),intent(in)    :: frac_npp_2_exudates    ! Fraction of NPP to be put into the optimally into the root exudates
    real(wp),intent(in)    :: frac_green_2_herbivory
    real(wp),intent(in)    :: tau_Cpool_litter_green ! Time constant by which green litter pool is depreciated  [days]
    real(wp),intent(in)    :: tau_Cpool_litter_wood  ! Time constant by which woody litter pool is depreciated  [days]
    real(wp),intent(in)    :: tau_Cpool_woods        ! Time constant by which woods Pool is depreciated  !!! [years] !!!
    real(wp),intent(in)    :: LAI_shed_constant      ! Leaf shedding at a constant rate for evergreens etc. [1/day]
    real(wp),intent(in)    :: frac_C_litter_green2atmos ! Fraction of heterotrophic loss of green litter pools emitted to
                                                        ! the atmosphere
    real(wp),intent(in)    :: Max_C_content_woods    ! Maximum carbon content of wood pool (from lctLib)
    real(wp),intent(in)    :: specific_leaf_area_C   ! Specific leaf area (from lctLib) [m^2(leaf)/mol(Carbon)]
    real(wp),intent(in)    :: reserveC2leafC         ! Ratio of max. carbon in reserve pool to max. carbon of leaves
    !logical, intent(in)    :: is_vegetation          ! logical vegetation mask
    integer,intent(in)     :: pheno_type             ! crop =5
    integer,intent(in)     :: time_steps_per_day
    real(wp),intent(in)    :: dtime
    !
    real(wp),intent(inout) :: Cpool_green            ! Green carbon pool: on input last value; updated on output 
                                                     !    [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: Cpool_woods            ! Wood carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: Cpool_reserve          ! Reserve carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: Cpool_crop_harvest     ! Crop harvest carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    !
    real(wp),intent(out)   :: soilResp_rate          ! Soil (=heterotrophic) respiration rate  [mol(C)/m^2 s] Note: this is a loss,
                                                     !    therefore negative!
    !real(wp),intent(out)   :: soilResp_rate_pot      ! soil respiration without N limitation     
    real(wp),intent(out)   :: NPP_flux_correction    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                                     !    correction arises either because otherwise the reserve pool would get
                                                     !    negative (positive correction), or the wood pool would exceed its
                                                     !    maximum value (negative correction). [mol(C)/m^2 s]
    real(wp),intent(out)   :: excess_NPP             ! Part of NPP that could not be stored in one of the plant carbon pools 
                                                     !    (green, wood, reserve), but had to be  thrown away (a posteriori 
                                                     !    reduction of NPP) [mol(C)/m^2 s]
    real(wp),intent(out)   :: root_exudates          ! Value of NPP_2_rootExudates had to be dropped into the litter_green_bg pool.
                                                     !    [mol(C)/m^2 s]
    real(wp),intent(out)   :: Cflx_litterTotal       ! Total carbon flux from vegetation to litter (green+woody, ag+bg)
    real(wp),intent(out)   :: Cflx_herbivory         !
    real(wp),intent(out)   :: Cflx_herbivory_LG      !
    real(wp),intent(out)   :: Cflx_herbivory_2_atm   !
    real(wp),intent(out)   :: NPP_act_yDayMean       ! Actual NPP after N-limitation and excess carbon drop [mol(C)/m^2 s]
    !
    ! added for thonicke fire algorithm (spitfire) needed to correctly split the carbon pools into the different fuel classes
    !real(wp),optional,intent(inout) :: frac_litter_wood_new   ! new fraction in above ground wood litter pool  

    ! Meteorology for Yasso
    real(wp),optional,intent(in)    :: temp2_30d              ! 30 day mean temperature
    real(wp),optional,intent(in)    :: precip_30d             ! 30 day mean precipitation

    ! Yasso pools
    ! Size class 1; green litter
    !   Aboveground  
    real(wp),optional,intent(inout) :: YCpool_acid_ag1        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_water_ag1       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_ethanol_ag1     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_nonsoluble_ag1  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground  
    real(wp),optional,intent(inout) :: YCpool_acid_bg1        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_water_bg1       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_ethanol_bg1     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_nonsoluble_bg1  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_humus_1         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2; woody litter
    !   Aboveground  
    real(wp),optional,intent(inout) :: YCpool_acid_ag2        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_water_ag2       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_ethanol_ag2     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_nonsoluble_ag2  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground  
    real(wp),optional,intent(inout) :: YCpool_acid_bg2        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_water_bg2       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_ethanol_bg2     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_nonsoluble_bg2  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),optional,intent(inout) :: YCpool_humus_2         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Parameters for Yasso
       ! Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    real(wp),optional,intent(in)    :: LeafLit_coef1          ! fraction going to the acid soluble pools 
    real(wp),optional,intent(in)    :: LeafLit_coef2          ! fraction going to the water soluble pools
    real(wp),optional,intent(in)    :: LeafLit_coef3          ! fraction going to the ethanol soluble pools
    real(wp),optional,intent(in)    :: LeafLit_coef4          ! fraction going to the non soluble pools
    real(wp),optional,intent(in)    :: LeafLit_coef5          ! fraction going to the humus pool
       ! Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    real(wp),optional,intent(in)    :: WoodLit_coef1          ! fraction going to the acid soluble pools 
    real(wp),optional,intent(in)    :: WoodLit_coef2          ! fraction going to the water soluble pools
    real(wp),optional,intent(in)    :: WoodLit_coef3          ! fraction going to the ethanol soluble pools
    real(wp),optional,intent(in)    :: WoodLit_coef4          ! fraction going to the non soluble pools
    real(wp),optional,intent(in)    :: WoodLit_coef5          ! fraction going to the humus pool
    real(wp),optional,intent(in)    :: WoodLitterSize         ! size of coarse debris



    ! locals
    logical  :: is_crop                 ! logical crop mask

    real(wp) :: d_Cpool_slow
    real(wp) :: d_Cpool_litter_green_ag, d_Cpool_litter_green_bg
    real(wp) :: d_Cpool_litter_wood_ag, d_Cpool_litter_wood_bg
    real(wp) :: Cpool_green_max,Cpool_reserve_optimal
    real(wp) :: NPP_2_greenPool, NPP_2_woodPool, NPP_2_reservePool, NPP_2_rootExudates, Green_2_herbivory
    real(wp) :: C_2_litter_greenPool, C_2_litter_woodPool
    real(wp) :: excess_carbon
    real(wp) :: leaf_shedding
    real(wp) :: alpha_mod

    real(wp) :: Cpool_litter_green_ag  ! Above ground green litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(wp) :: Cpool_litter_green_bg  ! Below ground green litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(wp) :: Cpool_litter_wood_ag   ! Above ground woody litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(wp) :: Cpool_litter_wood_bg   ! Below ground woody litter C-pool: updated on output [mol(C)/m^2(canopy)]
    real(wp) :: Cpool_slow             ! Slow soil carbon pool: on input last value; updated on output 
                                                     ! [mol(C)/m^2(canopy)]

    real(wp) :: Cpool_green_pot, Cpool_woods_pot         !! values of potential carbon pools i.e. before accounting for N-limitation
    real(wp) :: Cpool_slow_pot, Cpool_reserve_pot
    real(wp) :: Cpool_litter_green_ag_pot, Cpool_litter_green_bg_pot
    real(wp) :: Cpool_litter_wood_ag_pot, Cpool_litter_wood_bg_pot

    real(wp) :: Cflx_NPP_2_green_pot,Cflx_NPP_2_wood_pot,Cflx_NPP_2_reserve_pot 
    real(wp) :: Cflx_litterGreen_ag_2_atmos_pot, Cflx_litterGreen_bg_2_atmos_pot
    real(wp) :: Cflx_litterGreen_ag_2_slow_pot, Cflx_litterGreen_bg_2_slow_pot
    real(wp) :: Cflx_litter_wood_ag_2_atmos_pot,Cflx_litter_wood_bg_2_atmos_pot
    real(wp) :: Cflx_litter_wood_ag_2_slow_pot,Cflx_litter_wood_bg_2_slow_pot
    real(wp) :: Cflx_2_crop_harvest, Cflx_crop_harvest_2_atmos
    real(wp) :: Cflx_slow_2_atmos
    real(wp) :: Cflx_faeces_2_LG, Cflx_faeces_2_atm
    ! R: only for spitfire
    !real(wp) :: litter_wood_new, litter_wood_total

    ! local variables for Yasso
    real(wp), dimension(18) :: Yasso_io_pools
    real(wp), dimension(2)  :: Weather
    real(wp), dimension(16) :: Yasso_out 
    real(wp), dimension(5)  :: LeafLit_coefV
    real(wp), dimension(5)  :: WoodLit_coefV   
    real(wp)                :: leafLitter      ! amount of leafLitter 
    real(wp)                :: sizeZero

    !real(wp) :: Nflux_NPP_2_green_pot,Nflux_NPP_2_wood_pot 
    !real(wp) :: Nflux_litterGreen_ag_2_slow_pot, Nflux_litterGreen_bg_2_slow_pot
    !real(wp) :: Nflux_litter_wood_ag_2_slow_pot,Nflux_litter_wood_bg_2_slow_pot
    !real(wp) :: minNflux_litter_green_ag_pot,minNflux_litter_green_bg_pot
    !real(wp) :: minNflux_litter_wood_ag_pot,minNflux_litter_wood_bg_pot
    !real(wp) :: N2O_emissions_litter_pot 

    !real(wp) :: minN_plant_demand,minN_soil_demand,total_minN_demand
    real(wp) :: decompRate_litter_green

    real(wp) :: NPP_2_green_via_mobileN, NPP_2_wood_via_mobileN
    !real(wp) :: Nmobile_2_green, Nmobile_2_wood, Nmobile_2_green_pot, Nmobile_2_wood_pot
    !real(wp) :: redFact_Nmobile
    real(wp) :: redFact_Nlimitation
    !real(wp) :: Nflx_faeces_2_SMINN, minNflux_litter_active
    !real(wp) :: Nflx_crop_harvest_2_SMINN

    !LOGICAL  :: with_Nitrogen        !! .FALSE.: only Carbon pools are updated; .TRUE.: in addition also Nitrogen pools are updated 
    !LOGICAL  :: with_yasso           !! .FALSE.: cbalance model for litter and soil carbon is used; .TRUE.: yasso model is used 
    !REAL(wp), PARAMETER :: N2O_ef_grazing = 0._wp  !! N2O grazing switched off; reasonable value for the emission factor: 0.0125_wp 


    ! R: Normally parameters should be in mo_carbon_constants, but the following parameters are used only in this subroutine!
    !    Therefore I put them here. Is there any not structural advantage to put them into mo_carbon_constants?
    ! Parameter for Subroutine calc_Cpools 
    REAL(wp), parameter :: Q10                    = 1.8_wp          ! Empirical parameter for temperature dependence of 
                                                                    !   heterotrophic respiration rate
    REAL(wp), parameter :: kappa                  = 1.0_wp          ! Empirical parameter for soil moisture dependence of 
                                                                    !   heterotrophic respiration rate
    REAL(wp), parameter :: referenceTemp_Q10      = Tmelt           ! Reference temperature in the Q10 formula [Kelvin]
    REAL(wp), parameter :: tau_Cpool_reserve      = 1.0_wp *365._wp ! time constant by which reserve pool is depreciated [days]

    REAL(wp), parameter :: tau_Cpool_slow         = 100._wp*365._wp ! time constant by which slow pool is depreciated  [days]
    REAL(wp), parameter :: tau_Cpool_crop_harvest = 1._wp*365._wp   ! time constant by which the crop harvest pool decays [days]
    REAL(wp), parameter :: greenC2leafC           = 4.0_wp          ! ratio of carbon green pool (leaves, fine roots, starches,
                                                                    !    sugars) to leaf carbon
    REAL(wp), parameter :: frac_C_litter_wood2atmos = 0.2_wp        ! fraction of heterotrophic loss of woody litter emitted to the
                                                                    !    atmosphere
    REAL(wp), parameter :: frac_C_faeces2_LG      = 0.3_wp          ! fraction of C from faces put into the litter green pool
    REAL(wp), parameter :: frac_C_crop_harvest    = 0.5_wp          ! fraction of C from ag litter flux put to the crop harvest pool
    REAL(wp), parameter :: alpha_critical         = 0.35_wp         ! critical value of relative soil moisture below which
                                                                    !    heterotrophic respiration almost stops
    REAL(wp), parameter :: alpha_min              = 0.1_wp          ! effective lowest value for soil moisture entering the
                                                                    !    heterotrophic respiration
    REAL(wp), parameter :: N2O_rate_nitrification = 0.001_wp        ! JSBACH assumes lower bound of estimation in the literture as: 
                                                                    !   (0.1%) due to nitrification 
    REAL(wp), parameter :: N2O_rate_denitrification = 0.00125_wp    ! 0.0025_wp  ! & (0.2%) due to denitrification 
                                                                    ! actually 0.2% of denitrified N, not of N inputs
                                                                    ! assumption: 50% loss of NO3 by leaching
    REAL(wp), parameter :: N2O_rate               = 0.005_wp        ! JSBACH assumes 0.5-1%
                                                                    !    N2O emissions rate per day (<0.1%-0.2% nitrification ); 
                                                                    !    0.2-4.7% due to denit in Xu Ri, 2008 * therein
    REAL(wp), parameter :: sminn_NH4_fraction     = 0.4_wp          ! soil mineral NH4 is 40% of SMINN pool and rest in the form of
                                                                    !  NO3. Ref: Xu Ri, 2008



   ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    soilResp_rate            = 0.0_wp
    !soilResp_rate_pot        = 0.0_wp
    NPP_act_yDayMean         = 0.0_wp
   ! redFact_Nlimitation      = 1.0_wp

    ! R: only for spitfire:
    !IF (PRESENT(frac_litter_wood_new)) THEN
    !   litter_wood_new = 0.0_wp
    !END IF

    ! Preparations
    !-------------------------------------------

    if (pheno_type == 5) THEN
       is_crop = .TRUE.
    ELSE
       is_crop = .FALSE.
    END IF

    ! modified soil moist. such that heterotrophic respiration is zero below alpha_critical
    alpha_mod = MAX(alpha_min, (alpha_yDayMean - alpha_critical) / (1._wp - alpha_critical))

    ! Initializations
    NPP_flux_correction   = 0.0_wp
    Cpool_green_max       = 0.0_wp
    Cpool_reserve_optimal = 0.0_wp
    excess_NPP            = 0.0_wp
    root_exudates         = 0.0_wp
    Cflx_litterTotal      = 0.0_wp
    Cflx_herbivory        = 0.0_wp
    Cflx_herbivory_LG     = 0.0_wp
    Cflx_herbivory_2_atm  = 0.0_wp
    Cflx_2_crop_harvest   = 0.0_wp

    ! with yasso:
       Yasso_io_pools(1)=YCpool_acid_ag1
       Yasso_io_pools(2)=YCpool_water_ag1
       Yasso_io_pools(3)=YCpool_ethanol_ag1
       Yasso_io_pools(4)=YCpool_nonsoluble_ag1
       Yasso_io_pools(5)=YCpool_acid_bg1
       Yasso_io_pools(6)=YCpool_water_bg1
       Yasso_io_pools(7)=YCpool_ethanol_bg1
       Yasso_io_pools(8)=YCpool_nonsoluble_bg1
       Yasso_io_pools(9)=YCpool_humus_1
       Yasso_io_pools(10)=YCpool_acid_ag2
       Yasso_io_pools(11)=YCpool_water_ag2
       Yasso_io_pools(12)=YCpool_ethanol_ag2
       Yasso_io_pools(13)=YCpool_nonsoluble_ag2
       Yasso_io_pools(14)=YCpool_acid_bg2
       Yasso_io_pools(15)=YCpool_water_bg2
       Yasso_io_pools(16)=YCpool_ethanol_bg2
       Yasso_io_pools(17)=YCpool_nonsoluble_bg2
       Yasso_io_pools(18)=YCpool_humus_2

       Weather(1) = temp2_30d
       Weather(2) = precip_30d

       LeafLit_coefV(1) = LeafLit_coef1
       LeafLit_coefV(2) = LeafLit_coef2
       LeafLit_coefV(3) = LeafLit_coef3
       LeafLit_coefV(4) = LeafLit_coef4
       LeafLit_coefV(5) = LeafLit_coef5
       WoodLit_coefV(1) = WoodLit_coef1
       WoodLit_coefV(2) = WoodLit_coef2
       WoodLit_coefV(3) = WoodLit_coef3
       WoodLit_coefV(4) = WoodLit_coef4
       WoodLit_coefV(5) = WoodLit_coef5

       Yasso_out(1:16) = 0.0_wp           ! array to store the Yasso pools/fluxes 
       leafLitter = 0.0_wp                ! amount of leafLitter 
       sizeZero = 0.0_wp


    ! on all veg tiles:

       Cpool_green_max = greenC2leafC * LAI / specific_leaf_area_C
       Cpool_reserve_optimal = reserveC2leafC * max_LAI / specific_leaf_area_C

       !! ============== NON-N-LIMITED CARBON AND NITROGEN ALLOCATION =========================================================
       !!
       !! Perform those C- and N-allocation steps that are not restricted by N-availabilty
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

       C_2_litter_greenPool =  &                                     !! Maximally the whole carbon from the green pool is transfered
          MIN(greenC2leafC * leaf_shedding / specific_leaf_area_C, & !!    by leaf and root shedding to the green litter pool
              Cpool_green)
       Cpool_green = Cpool_green - C_2_litter_greenPool              !! This removes carbon from the green pool
       excess_carbon = max(0.0_wp,Cpool_green - Cpool_green_max)     !! If by the reduction of LAI the maximum value of the green
                                                                     !!    pool is still smaller than the current value, there is 
       C_2_litter_greenPool = C_2_litter_greenPool + excess_carbon   !!    excess C that has also to be shedded to the green litter
       Cpool_green = Cpool_green - excess_carbon                     !!    pool and also substracted from the green pool

       IF (is_crop) THEN
          !!  Move shedded Carbon pro rata to crop harvest flux and to ag green litter pools
          Cflx_2_crop_harvest = frac_green_aboveGround * frac_C_crop_harvest * C_2_litter_greenPool
          Cpool_litter_green_ag = Cpool_litter_green_ag + &
             frac_green_aboveGround * (1._wp - frac_C_crop_harvest) * C_2_litter_greenPool
       ELSE
          !! Move shedded Carbon to ag green litter pools
          Cpool_litter_green_ag = Cpool_litter_green_ag + frac_green_aboveGround * C_2_litter_greenPool
       END IF
       !! Move shedded Carbon to bg green litter pools
       Cpool_litter_green_bg = Cpool_litter_green_bg + (1._wp - frac_green_aboveGround) * C_2_litter_greenPool


       Cflx_litterTotal = C_2_litter_greenPool

       !! ---------------------------------------------------------------------------------
       !! Herbivory loss from Cpool_green & Npool_green w.r.t PFTs in lctlib

       Green_2_herbivory  = Cpool_green * frac_green_2_herbivory             !! Grazing flux indirect comput.from NPP losses 
       Cpool_green        = Cpool_green - Green_2_herbivory                  !! Loss due to grazing

       Cflx_faeces_2_LG  = Green_2_herbivory * frac_C_faeces2_LG             !! Cflux_faeces_2_LG to be stored in litter green pool
                                                                             !!      (not considered in Cflx_litterTotal) 
                                                                             !!      (no life times (tau) enter Cflx_faeces_2_LG!)
       Cflx_faeces_2_atm = Green_2_herbivory * (1.0_wp - frac_C_faeces2_LG)  !! Cflux_faeces_2_atm to be lost to atmosphere 
                                                                             !!      (dk: not added to CO2?)
       Cpool_litter_green_ag = Cpool_litter_green_ag + Cflx_faeces_2_LG      !! Cflux_faeces_2_LG put into LG pool  



       !! ---------------------------------------------------------------------------------
       !! 2. Wood shedding: Transfer of C from wood pool to wood litter pools
       !! 
                                                                                     !! Assuming that forests continously die at the inverse lifetime of trees
       C_2_litter_woodPool = Cpool_woods /              &                            ! The bigger tau gets, the weaker the C pools will be depreciated:
                                                                                     ! tau_Cpool_woods: [years] -> [timestep]
                             &  (tau_Cpool_woods * 365._wp )                         !! .. the shedded wood (MAX FUNCTION FOR NONWOODY PFTS)
       Cpool_litter_wood_ag = Cpool_litter_wood_ag      &                            !! .. is put partly into the above ground woody litter pool
            + frac_wood_aboveGround*C_2_litter_woodPool
       ! R: only for spitfire
       ! if (present(frac_litter_wood_new))     litter_wood_new = frac_wood_aboveGround * C_2_litter_woodPool  
       Cpool_litter_wood_bg = Cpool_litter_wood_bg &                    !! .. and the rest into the below ground woody litter pool
            + (1.0_wp-frac_wood_aboveGround)*C_2_litter_woodPool
       Cpool_woods = Cpool_woods -  C_2_litter_woodPool                 !! .. and then all is substracted from the wood pool.

       Cflx_litterTotal = Cflx_litterTotal + C_2_litter_woodPool



       !! ----------------------------------------------------------------------------
       !! 3. Depletion of reserve pool: Transfer of C from reserve pool to green litter pool
       !!    Some organic carbon of the reserve pool is always lost with mortality of plants. 
       !!    Note that the reserve pool contains no Nitrogen (starches and sugar are free of N).
       !!    Note also that the green litter pool has no fixed C/N-ratio (therefore no compensation
       !!    flux is needed)

       C_2_litter_greenPool = Cpool_reserve / tau_Cpool_reserve
       IF (is_crop) THEN
          Cflx_2_crop_harvest = Cflx_2_crop_harvest + &
              frac_green_aboveGround * frac_C_crop_harvest * C_2_litter_greenPool
          Cpool_litter_green_ag = Cpool_litter_green_ag + &
              frac_green_aboveGround * (1._wp - frac_C_crop_harvest) * C_2_litter_greenPool
       ELSE
          Cpool_litter_green_ag = Cpool_litter_green_ag + frac_green_aboveGround * C_2_litter_greenPool
       END IF
       Cpool_litter_green_bg = Cpool_litter_green_bg + (1._wp - frac_green_aboveGround) * C_2_litter_greenPool
       Cpool_reserve = Cpool_reserve - C_2_litter_greenPool
       Cpool_crop_harvest = Cpool_crop_harvest + Cflx_2_crop_harvest

       Cflx_litterTotal = Cflx_litterTotal + C_2_litter_greenPool

       !! ----------------------------------------------------------------------------
       !! 4. Decomposition of slow soil pool (i) C-flux (ii) Mineral N-flux
       !!
       d_Cpool_slow = &
            -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._wp)  / & !! Slow soil respiration
            tau_Cpool_slow,1._wp)  * Cpool_slow                                          !! .. according to Q10-model
       Cpool_slow = Cpool_slow + d_Cpool_slow                 !! Remove respired Carbon (note the negative sign!)

       Cflx_slow_2_atmos            = -d_Cpool_slow          !! Remember the decomposition flux to the atmosphere (CO2-release)


       !!
       !! ============== START OF POTENTIAL PLANT-CARBON ALLOCATION TO ESTIMATE PLANT-N-DEMAND =================================
       !!                       (exception: case NPP<0 handles actual allocation)

       Cpool_green_pot   = Cpool_green
       Cpool_woods_pot   = Cpool_woods
       Cpool_reserve_pot = Cpool_reserve

       !! Preliminary determination of distribution of NPP to the various pools (may be corrected afterwards according to 
       !!    avaialable Nitrogen)

       if(NPP_pot_yDayMean <= 0.0_wp) then                             !! In case of negative or zero NPP
          NPP_2_reservePool =  NPP_pot_yDayMean * sec_per_day          !! .. the plant looses C and  
          Cpool_reserve     = Cpool_reserve + NPP_2_reservePool        !! .. we try to take it from the reserve pool 
          if(Cpool_reserve < 0.0_wp) then                              !! But if thereby the reserve pool gets negative 
             NPP_2_reservePool = NPP_2_reservePool - Cpool_reserve     !! .. we absorb negative NPP only up to zero reseve pool and 
             NPP_flux_correction = -Cpool_reserve / sec_per_day        !! .. give the deficit as a flux correction back to the
             Cpool_reserve = 0.0_wp                                    !! .. calling routine. Hence all Carbon from the reserve pool
          end if                                                       !! .. is used up.

          !! For correct handling of nitrogen limitation below (which handles only the case of positive NPP) set NPP to zero:
          NPP_2_reservePool = 0.0_wp
          NPP_2_greenPool   = 0.0_wp   !! All other transfer rates are zero
          NPP_2_woodPool    = 0.0_wp
          NPP_2_rootExudates= 0.0_wp
          Cpool_reserve_pot = Cpool_reserve      !! remember modified reserve pool


       else  !! NPP is positive                                                            

          NPP_2_woodPool    = frac_npp_2_woodPool * NPP_pot_yDayMean * sec_per_day  !! NPP it is distributed according to predefined 
                                                                                    !!    relative fraction
          NPP_2_reservePool = frac_npp_2_reservePool * NPP_pot_yDayMean * sec_per_day
          NPP_2_rootExudates= frac_npp_2_exudates * NPP_pot_yDayMean * sec_per_day
          NPP_2_greenPool   = (1.0_wp - frac_npp_2_woodPool - frac_npp_2_reservePool                 & 
               - frac_npp_2_exudates) * NPP_pot_yDayMean * sec_per_day

          !! Growth of Wood Pool (structural carbon of living plants)

          Cpool_woods_pot = Cpool_woods_pot + NPP_2_woodPool                 !! Then it is attempted to put all NPP share into it
          excess_carbon = max(0.0_wp,Cpool_woods_pot-Max_C_content_woods)    !! .. but thereby the pool may get too large
          Cpool_woods_pot = Cpool_woods_pot - excess_carbon                  !! .. so that the excess carbon has once more to be 
                                                                             !! .. subtracted
          NPP_2_greenPool = NPP_2_greenPool + excess_carbon                  !! .. and instead made available to the green pool
          NPP_2_woodPool = NPP_2_woodPool - excess_carbon                    !! .. and the actual amount of NPP put to the wood pool
                                                                             !! .. is less


          !! Growth of Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
          if(Cpool_reserve_pot < Cpool_reserve_optimal) then                 !! .. If the reserve pool is smaller than optimal,
             Cpool_reserve_pot = Cpool_reserve_pot +  NPP_2_reservePool      !! .... it is filled by the available NPP
             excess_carbon =                                         &       !! .... Thereby it may happen that it gets larger than
                  max(0.0_wp,Cpool_reserve_pot - Cpool_reserve_optimal)      !! .... optimal so that there is excess carbon
             Cpool_reserve_pot = Cpool_reserve_pot - excess_carbon           !! .... that needs not be taken up
             NPP_2_greenPool = NPP_2_greenPool + excess_carbon               !! .... but can better be used to increase the green
             NPP_2_reservePool = NPP_2_reservePool - excess_carbon           !! .... pool and the actual amount of NPP put to the
                                                                             !! .... reserve pool is less.
          else                                                               !! .... Otherwise (reserve pool is larger than optimal)
             NPP_2_greenPool = NPP_2_greenPool + NPP_2_reservePool           !! .... all NPP is left for the green pool
             NPP_2_reservePool = 0.0_wp                                      !! .... so that nothing is stored in the reserve pool.
          end if

          !! Growth of Green Pool (leaves and fine roots): (in case of too much NPP, try to put it into the reserve pool).

          Cpool_green_pot = Cpool_green_pot + NPP_2_greenPool                !! Green pool is filled by the available NPP.
          excess_carbon = max(0.0_wp,Cpool_green_pot - Cpool_green_max)      !! .. Thereby it may get larger than appropriate for 
                                                                             !!    current LAI.
          NPP_2_greenPool = NPP_2_greenPool - excess_carbon                  !! .. Hence the actual amount of NPP put to the green
                                                                             !!    pool is less.
          Cpool_green_pot = Cpool_green_pot - excess_carbon                  !! .. This excess carbon needs not be taken up but
          if(Cpool_reserve_pot < Cpool_reserve_optimal) then                 !! .... if the reserve pool is smaller than optimal
             Cpool_reserve_pot = Cpool_reserve_pot + excess_carbon           !! .... it is tried to put the carbon there,
             NPP_2_reservePool = NPP_2_reservePool + excess_carbon           !! .... which means that additional NPP is put to the 
                                                                             !!      reserve pool.
             excess_carbon =                                          &      !! .... Thereby it may happen that the reserve pool
                  max(0.0_wp,Cpool_reserve_pot - Cpool_reserve_optimal)      !!      increases beyond the optimal value.
             Cpool_reserve_pot = Cpool_reserve_pot - excess_carbon           !! .... In that case the excess carbon is once more
                                                                             !!      removed from the reserve pool
             NPP_2_reservePool = NPP_2_reservePool - excess_carbon           !! .... so that the actual amount of NPP put into the
                                                                             !!      reserve pool is less.
          end if                                                             !!

          Cpool_litter_green_bg = Cpool_litter_green_bg + NPP_2_rootExudates
          ! THT 30.10.2012: This taken into account in Yasso too.. Put NPP_2_rootExudates to Yasso-water-pool. Added 20.11.2012
          ! Yasso_io_pools(2) = Yasso_io_pools(2) + NPP_2_rootExudates ! DSG: exudates goes into yasso

          excess_NPP = excess_carbon / sec_per_day

          root_exudates = NPP_2_rootExudates / sec_per_day

       end if !! NPP > 0 end


       !! ============== START OF POTENTIAL LITTER-CARBON ALLOCATION ==============================================

       Cpool_litter_green_ag_pot   = Cpool_litter_green_ag
       Cpool_litter_green_bg_pot   = Cpool_litter_green_bg
       Cpool_litter_wood_ag_pot    = Cpool_litter_wood_ag
       Cpool_litter_wood_bg_pot    = Cpool_litter_wood_bg

       ! Update Green Litter Pools (litter from dead leaves, fruits, debris from bark and fine roots)

       decompRate_litter_green =  &                                     !! decomposition rate of green litter according to Q10 model
            -MIN(alpha_mod**kappa * q10**((topSoilTemp - referenceTemp_Q10)/10._wp)/ tau_Cpool_litter_green,1._wp)

       d_Cpool_litter_green_ag      = decompRate_litter_green * Cpool_litter_green_ag_pot     !! Decomposition of green plant litter
       Cpool_litter_green_ag_pot    = Cpool_litter_green_ag_pot + d_Cpool_litter_green_ag
       d_Cpool_litter_green_bg      = decompRate_litter_green * Cpool_litter_green_bg_pot     !! Decomposition of green plant litter
       Cpool_litter_green_bg_pot    = Cpool_litter_green_bg_pot + d_Cpool_litter_green_bg

       ! Update Wood Litter Pools (litter from all dead woody material, above and below ground)

       d_Cpool_litter_wood_ag = -MIN(q10**((topSoilTemp - referenceTemp_Q10)/10._wp)  / &   !! Decomposition of woody litter
            tau_Cpool_litter_wood,1._wp) * Cpool_litter_wood_ag_pot                         !! .. according to Q10-model
       d_Cpool_litter_wood_bg = -MIN(q10**((topSoilTemp - referenceTemp_Q10)/10._wp)  / &   !! Decomposition of woody litter
            tau_Cpool_litter_wood,1._wp) * Cpool_litter_wood_bg_pot                         !! .. according to Q10-model
       !! .. without soil moisture dependence
       Cpool_litter_wood_ag_pot = Cpool_litter_wood_ag_pot + d_Cpool_litter_wood_ag
       Cpool_litter_wood_bg_pot = Cpool_litter_wood_bg_pot + d_Cpool_litter_wood_bg


       !! ============== START OF POTENTIAL ALLOCATION OF SLOW SOIL POOL ==============================================


       ! Update Slow Soil Pool (uptake of C from litter pools, decomposition has already been handled above)
       Cpool_slow_pot = Cpool_slow  
       Cpool_slow_pot = Cpool_slow_pot                                                           &
            - (1._wp - frac_C_litter_wood2atmos) * (d_Cpool_litter_wood_ag + d_Cpool_litter_wood_bg)  &
            - (1._wp - frac_C_litter_green2atmos) * (d_Cpool_litter_green_ag + d_Cpool_litter_green_bg)           



       !! ============== use nitrogen from N mobile pool to satisfy (at least partly) for potential carbon fluxes ===
       ! Remarks: Nmobile_2_green + Nmobile_2_wood <= Npool_mobile
       ! Npool_mobile accounts minN_plant_demand by trasfering N into green and wood pool 

       Cflx_NPP_2_green_pot          = NPP_2_greenPool                                       !! potential growth of green carbon 
       Cflx_NPP_2_wood_pot           = NPP_2_woodPool                                        !! potential growth of wood carbon 

       NPP_2_green_via_mobileN = 0.0_wp
       NPP_2_wood_via_mobileN = 0.0_wp 



       Cflx_NPP_2_green_pot = Cflx_NPP_2_green_pot - NPP_2_green_via_mobileN
       Cflx_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot  - NPP_2_wood_via_mobileN

       ! R: These lines were calculated also with_Nitrogen=.FALSE. in JS3. However, it was not used further when with_Nitrogen=.FALSE.
       !    In JS3 it should stand within the brackets.
       !Nflux_NPP_2_green_pot = Cflx_NPP_2_green_pot/cn_green
       !Nflux_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot/cn_woods
       !minN_plant_demand = Nflux_NPP_2_green_pot + Nflux_NPP_2_wood_pot 


       !! ============== COLLECT POTENTIAL C-FLUXES ===============================================================

       Cflx_NPP_2_reserve_pot        = NPP_2_reservePool                                     !! potential growth of reserve carbon
       Cflx_litterGreen_ag_2_atmos_pot = -frac_C_litter_green2atmos * d_Cpool_litter_green_ag !! ag green litter decomp. -> atm.
       Cflx_litterGreen_bg_2_atmos_pot = -frac_C_litter_green2atmos * d_Cpool_litter_green_bg !! bg green litter decomp. -> atm.
       Cflx_litterGreen_ag_2_slow_pot  = &
           - (1._wp -frac_C_litter_green2atmos) * d_Cpool_litter_green_ag                    !! ag green litter decomp. -> slow pool
       Cflx_litterGreen_bg_2_slow_pot  = &
            -(1._wp -frac_C_litter_green2atmos) * d_Cpool_litter_green_bg                    !! bg green litter decomp. -> slow pool
       Cflx_litter_wood_ag_2_atmos_pot  = -frac_C_litter_wood2atmos * d_Cpool_litter_wood_ag !! ag-wood litter decomp.  -> atm.
       Cflx_litter_wood_bg_2_atmos_pot  = -frac_C_litter_wood2atmos * d_Cpool_litter_wood_bg !! bg-wood litter decomp.  -> atm.
       Cflx_litter_wood_ag_2_slow_pot   = &
            -(1._wp - frac_C_litter_wood2atmos) * d_Cpool_litter_wood_ag                     !! ag-wood litter decomp.  -> slow pool
       Cflx_litter_wood_bg_2_slow_pot   = &
            -(1._wp - frac_C_litter_wood2atmos) * d_Cpool_litter_wood_bg                     !! bg_wood litter decomp.  -> slow pool

 



          redFact_Nlimitation = 1.0_wp   !! This assures absence of nitrogen limitation in this mixed C-N code


       !! ============== UPDATE ALL C-POOLS =======================================================================

       Cpool_green       = Cpool_green       + redFact_Nlimitation * Cflx_NPP_2_green_pot
       Cpool_woods       = Cpool_woods       + redFact_Nlimitation * Cflx_NPP_2_wood_pot
       Cpool_reserve     = Cpool_reserve     +  Cflx_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N

       Cpool_litter_green_ag = Cpool_litter_green_ag       &
            - redFact_Nlimitation * (Cflx_litterGreen_ag_2_atmos_pot + Cflx_litterGreen_ag_2_slow_pot)

       Cpool_litter_green_bg = Cpool_litter_green_bg       &
            - redFact_Nlimitation * (Cflx_litterGreen_bg_2_atmos_pot + Cflx_litterGreen_bg_2_slow_pot)

       Cpool_litter_wood_ag = Cpool_litter_wood_ag       &
            - redFact_Nlimitation * (Cflx_litter_wood_ag_2_atmos_pot + Cflx_litter_wood_ag_2_slow_pot)

       Cpool_litter_wood_bg = Cpool_litter_wood_bg       &
            - redFact_Nlimitation * (Cflx_litter_wood_bg_2_atmos_pot + Cflx_litter_wood_bg_2_slow_pot)

       Cpool_slow        = Cpool_slow     &
            + redFact_Nlimitation * (  Cflx_litterGreen_ag_2_slow_pot + Cflx_litterGreen_bg_2_slow_pot &
                                     + Cflx_litter_wood_ag_2_slow_pot  + Cflx_litter_wood_bg_2_slow_pot  )





       !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

       if(NPP_pot_yDayMean > 0.0_wp) then 
          NPP_act_yDayMean = (NPP_2_green_via_mobileN + NPP_2_wood_via_mobileN + NPP_2_rootExudates +               &
               redFact_Nlimitation * (Cflx_NPP_2_green_pot + Cflx_NPP_2_wood_pot) + Cflx_NPP_2_reserve_pot &
               )/sec_per_day
       else
          NPP_act_yDayMean = NPP_pot_yDayMean + NPP_flux_correction  ! NPP_pot_yDayMean is negative: NPP_flux_correction is zero 
                                                                     ! except Cpool_reserve was < 0.0_wp.
       endif


       !! ============== COMPUTE SOIL RESPIRATION ===================================================================
       ! with yasso:
          leafLitter    = Cflx_litterTotal - C_2_litter_woodPool + Cflx_faeces_2_LG - Cflx_2_crop_harvest
          soilResp_rate = 0._wp

          !! call yasso for leaf litter; this time with nutrient effects
          CALL yasso(Yasso_io_pools(1:9), Weather, leafLitter, LeafLit_coefV,               & 
                           sizeZero, Yasso_out, frac_green_aboveGround, NPP_2_rootExudates,redFact_Nlimitation)

          YCpool_acid_ag1         =Yasso_out(1)
          YCpool_water_ag1        =Yasso_out(2)
          YCpool_ethanol_ag1      =Yasso_out(3)
          YCpool_nonsoluble_ag1   =Yasso_out(4)
          YCpool_acid_bg1         =Yasso_out(5)
          YCpool_water_bg1        =Yasso_out(6)
          YCpool_ethanol_bg1      =Yasso_out(7)
          YCpool_nonsoluble_bg1   =Yasso_out(8)
          YCpool_humus_1          =Yasso_out(9)


          ! rescale respiration (from leaf litter decomposition) from "per day" to "per second"
          soilResp_rate = Yasso_out(10) / sec_per_day 

          !! call yasso for woody litter; this time with nutrient effects
          CALL yasso(Yasso_io_pools(10:18), Weather, C_2_litter_woodPool, WoodLit_coefV,       & 
                           WoodLitterSize, Yasso_out, frac_wood_aboveGround, 0.0_wp, redFact_Nlimitation) ! no exudates
  
          YCpool_acid_ag2         =Yasso_out(1)
          YCpool_water_ag2        =Yasso_out(2)
          YCpool_ethanol_ag2      =Yasso_out(3)
          YCpool_nonsoluble_ag2   =Yasso_out(4)
          YCpool_acid_bg2         =Yasso_out(5)
          YCpool_water_bg2        =Yasso_out(6)
          YCpool_ethanol_bg2      =Yasso_out(7)
          YCpool_nonsoluble_bg2   =Yasso_out(8)
          YCpool_humus_2          =Yasso_out(9)
       
          !! Yasso total respiration
          ! add respiration from wood litter decomposition and rescale respiration from "per day" to "per second"
           soilResp_rate = soilResp_rate  + Yasso_out(10) / sec_per_day

       ! always:
       !! ============== Decay of crop harvest =============

       Cflx_crop_harvest_2_atmos = Cpool_crop_harvest / tau_Cpool_crop_harvest
       Cpool_crop_harvest = Cpool_crop_harvest - Cflx_crop_harvest_2_atmos
       soilResp_rate = soilResp_rate - Cflx_crop_harvest_2_atmos / sec_per_day
       !soilResp_rate_pot = soilResp_rate_pot - Cflx_crop_harvest_2_atmos / sec_per_day


       !! ============== COMPUTE Herbivory RESPIRATION =====

       Cflx_herbivory_2_atm = -(Cflx_faeces_2_atm) /sec_per_day

       !! ============== Litter and herbivory diagnostics (flux from day-1 to s-1)
       Cflx_litterTotal = Cflx_litterTotal /sec_per_day
       Cflx_herbivory    = Green_2_herbivory / sec_per_day
       Cflx_herbivory_LG = Cflx_faeces_2_LG / sec_per_day   ! DSG: goes into yasso



! R: only for spitfire:
       ! fraction of new wood litter, assume old and new litter are respired in the same way, respiration does not change the fraction
!        if (present(frac_litter_wood_new)) then
!           if (with_yasso) then
!             litter_wood_total = YCpool_water_ag2 + YCpool_acid_ag2 + YCpool_ethanol_ag2 + YCpool_nonsoluble_ag2
!           else
!             litter_wood_total = Cpool_litter_wood_ag
!           endif
!           if ((litter_wood_total >  EPSILON(1._wp)) .AND. (litter_wood_new < litter_wood_total)) then
!              frac_litter_wood_new = litter_wood_new / litter_wood_total
!           else
!              frac_litter_wood_new = 1._wp    ! frac_litter_wood_new=1 means fuel fractions will be reset to original values
!           endif
!        endif


  end subroutine calc_Cpools




  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of fire for the carbon pools. More precisely:
  ! It is assumed that for the burned area the carbon from the above ground litter pools (Cpool_litter_green_ag,
  ! Cpool_litter_wood_ag) is completely released to the atmosphere and the carbon from the living plant pools
  ! (Cpool_green_st, Cpool_reserve_st, Cpool_woods_st) is partly released to the atmosphere and partly relocated
  ! into the litter pools.
  !
  subroutine relocate_carbon_fire( &
    & burned_frac,                  & ! in from distrubed_frac
    & fire_frac_wood_2_atmos,       & ! in
    & LeafLit_coef,                 & ! in
    & WoodLit_coef,                 & ! in
    & Cpool_green,                  & ! inout
    & Cpool_woods,                  & ! inout
    & Cpool_reserve,                & ! inout
    & carbon_2_GreenLitterPools,    & ! inout
    & carbon_2_WoodLitterPools,     & ! inout
    & YCpool_acid_ag1,              & ! inout
    & YCpool_water_ag1,             & ! inout
    & YCpool_ethanol_ag1,           & ! inout
    & YCpool_nonsoluble_ag1,        & ! inout
    & YCpool_acid_ag2,              & ! inout
    & YCpool_water_ag2,             & ! inout
    & YCpool_ethanol_ag2,           & ! inout
    & YCpool_nonsoluble_ag2,        & ! inout
    & YCpool_acid_bg1,              & ! inout
    & YCpool_water_bg1,             & ! inout
    & YCpool_ethanol_bg1,           & ! inout
    & YCpool_nonsoluble_bg1,        & ! inout
    & YCpool_acid_bg2,              & ! inout
    & YCpool_water_bg2,             & ! inout
    & YCpool_ethanol_bg2,           & ! inout
    & YCpool_nonsoluble_bg2,        & ! inout
    & YCpool_humus_1,               & ! inout
    & YCpool_humus_2,               & ! inout
    & carbon_2_atmos                & ! out
    & )

    USE mo_carbon_constants,  ONLY : frac_green_aboveGround, frac_wood_aboveGround

    ! Arguments
    real(wp),intent(in)    :: burned_frac(:)              ! Fraction of the vegetated area of each tile burned till the last
                                                               ! call of this routine
    real(wp),intent(in)    :: fire_frac_wood_2_atmos           ! Fraction of above ground wood immediately emitted to atmosphere 
                                                               ! by fire
    real(wp),intent(in)    :: LeafLit_coef(:)           ! Factor to spread non woody litterfall into yasso pools [ ]
    real(wp),intent(in)    :: WoodLit_coef(:)           ! Factor to spread woody litterfall into yasso pools [ ]
    real(wp),intent(inout) :: Cpool_green(:)              ! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: Cpool_woods(:)              ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: Cpool_reserve(:)            ! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: carbon_2_GreenLitterPools(:) ! Amount of carbon relocated by wind and fire damage
                                                                ! .. to the green litter pools [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: carbon_2_WoodLitterPools(:)  ! Amount of carbon relocated by wind damage and fire
                                                                ! .. to the wood litter pools [mol(C)/m^2(canopy)]
    !   YASSO
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    real(wp),intent(inout) :: YCpool_acid_ag1(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_water_ag1(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_ethanol_ag1(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_nonsoluble_ag1(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: YCpool_acid_bg1(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_water_bg1(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_ethanol_bg1(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_nonsoluble_bg1(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_humus_1(:)           ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    real(wp),intent(inout) :: YCpool_acid_ag2(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_water_ag2(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_ethanol_ag2(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_nonsoluble_ag2(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(wp),intent(inout) :: YCpool_acid_bg2(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_water_bg2(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_ethanol_bg2(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_nonsoluble_bg2(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_humus_2(:)           ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    real(wp),intent(out)   :: carbon_2_atmos(:)           ! Carbon immedieatly released by fire per tile [mol(C)m-2(canopy)]

    ! Local variables
    REAL(wp) :: factor(SIZE(burned_frac)) ! Helper

    !! diagnose amount of carbon released from the green and reserve pool to the green litter pools
    !! diagnose amount of carbon released from the wood pool to the woody litter pools
    carbon_2_GreenLitterPools(:) = carbon_2_GreenLitterPools(:) + &
      & (Cpool_green(:) + Cpool_reserve(:)) * (1._wp - frac_green_aboveGround) * burned_frac(:)
    carbon_2_WoodLitterPools(:)  = carbon_2_WoodLitterPools(:) + &
      & Cpool_woods(:) * (1._wp - fire_frac_wood_2_atmos * frac_wood_aboveGround) * burned_frac(:)

    !! determine amount of carbon released from the wood pool, living tissue pools, and litter pools to the atmosphere
    carbon_2_atmos(:) = ( & 
        ! take living parts of plants above ground:
      & (Cpool_green(:) + Cpool_reserve(:)) * frac_green_aboveGround             &
        ! now add the YASSO pools above ground:
      & + YCpool_water_ag1(:) + YCpool_acid_ag1(:) + YCPool_ethanol_ag1(:)  &
      & + YCpool_nonsoluble_ag1(:)                                                    &
      & + YCpool_water_ag2(:) + YCpool_acid_ag2(:) + YCPool_ethanol_ag2(:)  &
      & + YCpool_nonsoluble_ag2(:)                                                    &
        ! now add the fire affected fraction of wood pools above ground:
      & + Cpool_woods(:) * fire_frac_wood_2_atmos * frac_wood_aboveGround             &
      & ) * burned_frac(:)


    ! Lower down density of aboveground litter pools
    factor(:) = 1._wp - burned_frac(:)
    YCpool_acid_ag1      (:) = factor(:) * YCpool_acid_ag1       (:)
    YCpool_water_ag1     (:) = factor(:) * YCpool_water_ag1      (:)
    YCpool_ethanol_ag1   (:) = factor(:) * YCpool_ethanol_ag1    (:)
    YCpool_nonsoluble_ag1(:) = factor(:) * YCpool_nonsoluble_ag1 (:)

    YCpool_acid_ag2      (:) = factor(:) * YCpool_acid_ag2       (:)
    YCpool_water_ag2     (:) = factor(:) * YCpool_water_ag2      (:)
    YCpool_ethanol_ag2   (:) = factor(:) * YCpool_ethanol_ag2    (:)
    YCpool_nonsoluble_ag2(:) = factor(:) * YCpool_nonsoluble_ag2 (:)

    ! Add green and reserve vegetation carbon to belowground green litter pools
    factor(:) = (1._wp-frac_green_aboveGround) * burned_frac(:) * (Cpool_green(:) + Cpool_reserve(:))
    YCpool_acid_bg1      (:) = YCpool_acid_bg1      (:) + factor(:) * LeafLit_coef(1)
    YCpool_water_bg1     (:) = YCpool_water_bg1     (:) + factor(:) * LeafLit_coef(2)
    YCpool_ethanol_bg1   (:) = YCpool_ethanol_bg1   (:) + factor(:) * LeafLit_coef(3)
    YCpool_nonsoluble_bg1(:) = YCpool_nonsoluble_bg1(:) + factor(:) * LeafLit_coef(4)
    YCpool_humus_1       (:) = YCpool_humus_1       (:) + factor(:) * LeafLit_coef(5)

    ! Add the remaings of burned wood to the aboveground wood pools:
    factor(:) = (1._wp-fire_frac_wood_2_atmos) * frac_wood_aboveGround * burned_frac(:) * Cpool_woods(:)
    YCpool_acid_ag2      (:) = YCpool_acid_ag2      (:) + factor(:) * WoodLit_coef(1)
    YCpool_water_ag2     (:) = YCpool_water_ag2     (:) + factor(:) * WoodLit_coef(2)
    YCpool_ethanol_ag2   (:) = YCpool_ethanol_ag2   (:) + factor(:) * WoodLit_coef(3) 
    YCpool_nonsoluble_ag2(:) = YCpool_nonsoluble_ag2(:) + factor(:) * WoodLit_coef(4) 
    YCpool_humus_2       (:) = YCpool_humus_2       (:) + factor(:) * WoodLit_coef(5)

    ! Add the remaing of burned wood to the belowground wood pools:
    factor(:) = (1._wp-frac_wood_aboveGround) * burned_frac(:) * Cpool_woods(:)
    YCpool_acid_bg2      (:) = YCpool_acid_bg2      (:) + factor(:) * WoodLit_coef(1)
    YCpool_water_bg2     (:) = YCpool_water_bg2     (:) + factor(:) * WoodLit_coef(2)
    YCpool_ethanol_bg2   (:) = YCpool_ethanol_bg2   (:) + factor(:) * WoodLit_coef(3)
    YCpool_nonsoluble_bg2(:) = YCpool_nonsoluble_bg2(:) + factor(:) * WoodLit_coef(4)
    YCpool_humus_2       (:) = YCpool_humus_2       (:) + factor(:) * WoodLit_coef(5)

    ! Lower down the carbon density of living plant pools
    factor(:) = 1._wp - burned_frac(:)
    Cpool_green(:)   = Cpool_green(:)   * factor(:)
    Cpool_woods(:)   = Cpool_woods(:)   * factor(:)
    Cpool_reserve(:) = Cpool_reserve(:) * factor(:)

  end subroutine relocate_carbon_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (Cpool_green_st, Cpool_reserve_st and 
  ! Cpool_woods_st) is partly relocated into the litter pools.
  !
  subroutine relocate_carbon_damage(                 &
    & damaged_frac,                                  & ! in
    & LeafLit_coef,                                  & ! in
    & WoodLit_coef,                                  & ! in
    & Cpool_green,                                   & ! inout, optional
    & Cpool_woods,                                   & ! inout, optional
    & Cpool_reserve,                                 & ! inout, optional
    & carbon_2_GreenLitterPools,                     & ! inout
    & carbon_2_WoodLitterPools,                      & ! inout
    & YCpool_acid_ag1,                               & ! inout
    & YCpool_water_ag1,                              & ! inout
    & YCpool_ethanol_ag1,                            & ! inout
    & YCpool_nonsoluble_ag1,                         & ! inout
    & YCpool_acid_ag2,                               & ! inout
    & YCpool_water_ag2,                              & ! inout
    & YCpool_ethanol_ag2,                            & ! inout
    & YCpool_nonsoluble_ag2,                         & ! inout
    & YCpool_acid_bg1,                               & ! inout
    & YCpool_water_bg1,                              & ! inout
    & YCpool_ethanol_bg1,                            & ! inout
    & YCpool_nonsoluble_bg1,                         & ! inout
    & YCpool_acid_bg2,                               & ! inout
    & YCpool_water_bg2,                              & ! inout
    & YCpool_ethanol_bg2,                            & ! inout
    & YCpool_nonsoluble_bg2,                         & ! inout
    & YCpool_humus_1,        YCpool_humus_2          & ! inout
    & )

    USE mo_carbon_constants,  ONLY : frac_green_aboveGround, frac_wood_aboveGround

    ! Arguments
    real(wp),intent(in)    :: damaged_frac(:)              ! Fraction of the vegetated area of each tile damaged till the last
                                                                ! call of this routine
    real(wp),intent(in)    :: LeafLit_coef(:)         ! Factor to spread non woody litter into yasso pools [ ]
    real(wp),intent(in)    :: WoodLit_coef(:)         ! Factor to spread woody litterfall into yasso pools [ ]
    real(wp),intent(inout) :: Cpool_green(:)               ! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: Cpool_woods(:)               ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: Cpool_reserve(:)             ! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(wp),intent(inout)   :: carbon_2_GreenLitterPools(:) ! Amount of carbon relocated by wind damage 
                                                                ! (and later also by fire)
    real(wp),intent(inout)   :: carbon_2_WoodLitterPools(:)  ! Amount of carbon relocated by wind damage 
                                                                ! (and later also by fire)
                                                                ! .. to the green litter pools [mol(C)/m^2(canopy)]

    !   YASSO
    ! Size class 1: green
    !   above ground C pools
    real(wp),intent(inout) :: YCpool_acid_ag1(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_water_ag1(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_ethanol_ag1(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: YCpool_nonsoluble_ag1(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(wp),intent(inout) :: YCpool_acid_bg1(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_water_bg1(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_ethanol_bg1(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: YCpool_nonsoluble_bg1(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_humus_1(:)         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2: wood 
    !   above ground C pools
    real(wp),intent(inout) :: YCpool_acid_ag2(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_water_ag2(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_ethanol_ag2(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: YCpool_nonsoluble_ag2(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools  
    real(wp),intent(inout) :: YCpool_acid_bg2(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_water_bg2(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_ethanol_bg2(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: YCpool_nonsoluble_bg2(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: YCpool_humus_2(:)         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp) :: factor(SIZE(damaged_frac))      ! Helper

    ! Initializations

    !! transfer carbon from the wood pool to the wood litter pools nd
    !! transfer carbon from the green and reserve pool to the green litter pool

    carbon_2_WoodLitterPools(:)  = carbon_2_WoodLitterPools(:) + &
      & Cpool_woods(:) * damaged_frac(:)
    carbon_2_GreenLitterPools(:) = carbon_2_GreenLitterPools(:) + &
      & (Cpool_green(:) + Cpool_reserve(:)) * damaged_frac(:)

    ! Add aboveground damaged green and reserve
    factor(:) = (Cpool_green(:) + Cpool_reserve(:)) * frac_green_aboveGround * damaged_frac(:)
    YCpool_acid_ag1      (:) = YCpool_acid_ag1      (:) + factor(:) * LeafLit_coef(1)
    YCpool_water_ag1     (:) = YCpool_water_ag1     (:) + factor(:) * LeafLit_coef(2)
    YCpool_ethanol_ag1   (:) = YCpool_ethanol_ag1   (:) + factor(:) * LeafLit_coef(3)
    YCpool_nonsoluble_ag1(:) = YCpool_nonsoluble_ag1(:) + factor(:) * LeafLit_coef(4)

    ! Add belowground damaged green and reserve
    factor(:) = (Cpool_green(:) + Cpool_reserve(:)) * (1._wp - frac_green_aboveGround) * damaged_frac(:)
    YCpool_acid_bg1      (:) = YCpool_acid_bg1      (:) + factor(:) * LeafLit_coef(1)
    YCpool_water_bg1     (:) = YCpool_water_bg1     (:) + factor(:) * LeafLit_coef(2)
    YCpool_ethanol_bg1   (:) = YCpool_ethanol_bg1   (:) + factor(:) * LeafLit_coef(3)
    YCpool_nonsoluble_bg1(:) = YCpool_nonsoluble_bg1(:) + factor(:) * LeafLit_coef(4)
    factor(:) = (Cpool_green(:) + Cpool_reserve(:)) * damaged_frac(:)
    YCpool_humus_1       (:) = YCpool_humus_1       (:) + factor(:) * LeafLit_coef(5)

    ! Add aboveground damaged wood
    factor(:) = frac_wood_aboveGround  * damaged_frac(:) * Cpool_woods(:)
    YCpool_acid_ag2      (:) = YCpool_acid_ag2      (:) + factor(:) * WoodLit_coef(1)
    YCpool_water_ag2     (:) = YCpool_water_ag2     (:) + factor(:) * WoodLit_coef(2)
    YCpool_ethanol_ag2   (:) = YCpool_ethanol_ag2   (:) + factor(:) * WoodLit_coef(3)
    YCpool_nonsoluble_ag2(:) = YCpool_nonsoluble_ag2(:) + factor(:) * WoodLit_coef(4)

    ! Add belowground damaged wood
    factor(:) = (1._wp - frac_wood_aboveGround)  * damaged_frac(:) * Cpool_woods(:)
    YCpool_acid_bg2      (:) = YCpool_acid_bg2      (:) + factor(:) * WoodLit_coef(1)
    YCpool_water_bg2     (:) = YCpool_water_bg2     (:) + factor(:) * WoodLit_coef(2)
    YCpool_ethanol_bg2   (:) = YCpool_ethanol_bg2   (:) + factor(:) * WoodLit_coef(3)
    YCpool_nonsoluble_bg2(:) = YCpool_nonsoluble_bg2(:) + factor(:) * WoodLit_coef(4)
    factor(:) = Cpool_woods(:) * damaged_frac(:)
    YCpool_humus_2       (:) = YCpool_humus_2       (:) + factor(:) * WoodLit_coef(5)

    ! Lower down the carbon density of living plant pools
    factor(:) = 1._wp - damaged_frac(:)
    Cpool_green(:)   = Cpool_green(:)   * factor(:)
    Cpool_woods(:)   = Cpool_woods(:)   * factor(:)
    Cpool_reserve(:) = Cpool_reserve(:) * factor(:)

  end subroutine relocate_carbon_damage



  PURE SUBROUTINE yasso (Yasso_io_pools, Weather, litter, Lit_coefV,WoodLitterSize, Yasso_out, frac_aboveground, &
                         NPP_2_rootExudates, redFact_Nlimit)
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
    !! distributed as part of the MPI Earth System Model under the MPI-M license conditions as stated in the header of this module. 
    !!
    !! model to calculate the amount of soil organic carbon, changes in the amount of soil organic carbon and
    !! heterotrophic soil respiration.
    !! For documention of yasso see Liski et al 2005, Tuomi et al 2008,2009,2011
    !!
    !! implementation: Tea Thum (FMI), Petri Räisänen (FMI), Daniel Goll (MPI-M)
    !!-------------------------------------------------------------------------------------------------------------------------

    USE mo_carbon_constants,       ONLY: days_per_year, sec_per_year
    USE mo_jsb_physical_constants, ONLY: tmelt
    IMPLICIT NONE

    REAL(wp), DIMENSION(9),  INTENT(IN)  :: Yasso_io_pools     ! Yasso pools IN                                        [mol(c)/m2]
    REAL(wp),                INTENT(IN)  :: WoodLitterSize     ! Size of woody litter; 
                                                               !     is zero when the subroutine is called for non-woody litter
    REAL(wp), DIMENSION(2),  INTENT(IN)  :: Weather            ! climatic drivers: air temperature and precipitation, 15-d means
                                                               ! Weather(1) in [K] ; Weather(2) in [mm/s]
    REAL(wp),                INTENT(IN)  :: litter             ! fresh litter inputs  (above and belowground)         [mol(C)/m2/day]
    REAL(wp), DIMENSION(5),  INTENT(IN)  :: Lit_coefV          ! fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(wp),                INTENT(IN)  :: frac_aboveground   ! parameter to seperate litter influx into above- and belowground 
                                                               !     part                                                   [ ]
    REAL(wp),                INTENT(IN)  :: NPP_2_rootExudates ! root exudates [daily]                                [mol(C)/m2/day]

    REAL(wp), DIMENSION(16), INTENT(OUT) :: Yasso_out          ! updated Yasso pools(1:9) [mol(c)/m2] & respiration(10) [mol(c)/m2/d]
                                                               ! & Cflx_2_humus (11) [mol(c)/m2/d] ... see update_Cpools for more 
                                                               ! information

    REAL(wp), OPTIONAL, INTENT(IN)       ::  redFact_Nlimit

    !local
    ! Yasso pools   [mol(c)/m2]
    REAL(wp)     :: YCpool_acid_ag    ! Aboveground  
    REAL(wp)     :: YCpool_water_ag  
    REAL(wp)     :: YCpool_ethanol_ag 
    REAL(wp)     :: YCpool_nonsoluble_ag
    REAL(wp)     :: YCpool_acid_bg    ! Belowground  
    REAL(wp)     :: YCpool_water_bg  
    REAL(wp)     :: YCpool_ethanol_bg 
    REAL(wp)     :: YCpool_nonsoluble_bg
    
    REAL(wp)     :: YCpool_humus   

    ! Yasso fluxes (only internally used and therefore as annual fluxes)  
    REAL(wp)     :: Cflx_litter_2_acid             ! Litter influx to acid soluble pools    [mol(c)/m2/a]   
    REAL(wp)     :: Cflx_litter_2_water            ! Litter influx to water soluble pools   [mol(c)/m2/a]   
    REAL(wp)     :: Cflx_litter_2_ethanol          ! Litter influx to ethanol soluble pools [mol(c)/m2/a]   
    REAL(wp)     :: Cflx_litter_2_nonsoluble       ! Litter influx to non soluble pools     [mol(c)/m2/a]   
    REAL(wp)     :: Cflx_litter_2_humus            ! Litter influx to humus soluble pool    [mol(c)/m2/a]   
  
    ! Aboveground
    REAL(wp)     :: Cflx_from_acid_ag              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_water_ag             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_ethanol_ag           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_nonsoluble_ag        ! Loss flux of non soluble pool         [mol(c)/m2/a]
    
    REAL(wp)     :: Cflx_2_acid_ag                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_water_ag                ! Gain flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_ethanol_ag              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_nonsoluble_ag           ! Gain flux of non soluble pool         [mol(c)/m2/a]
    
    ! Belowground
    REAL(wp)     :: Cflx_from_acid_bg              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_water_bg             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_ethanol_bg           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_nonsoluble_bg        ! Loss flux of non soluble pool         [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_humus                ! Loss flux of humus pool               [mol(c)/m2/a]
    
    REAL(wp)     :: Cflx_2_acid_bg                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_water_bg                ! Gain flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_ethanol_bg              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_nonsoluble_bg           ! Gain flux of non soluble pool         [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humus                   ! Gain flux of humus pool               [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humusAG                 ! Gain flux of humus pool; from ag only [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humusBG                 ! Gain flux of humus pool; from bg only [mol(c)/m2/a]
    
    ! Respiration is output of Yasso therefore a daily flux 
    REAL(wp)     :: soilResp_rateYasso             ! Flux to the atmosphere                [mol(c)/m2/d]  
    REAL(wp)     :: soilResp_rateLitterAG          ! Flux to the atmosphere; from AWEN ag only  [mol(c)/m2/d]
    REAL(wp)     :: soilResp_rateLitterBG          ! Flux to the atmosphere; from AWEN bg only  [mol(c)/m2/d]

    ! Decomposition rates [1/a]
    REAL(wp)     :: d_acid
    REAL(wp)     :: d_water
    REAL(wp)     :: d_ethanol
    REAL(wp)     :: d_nonsoluble
    REAL(wp)     :: d_humus

    ! Dcomposition rate for N litter green (needed for N cycle)
    REAL(wp)     :: d_litter_green              ! the rate with which the sum of AWEN pools decays
    REAL(wp)     :: Pseudo_litter_green         ! Sum of all AWEN pools

    ! PFT-specific parameters
    REAL(wp)     :: frac_litter_2_acid        ! Fraction of incoming material which enters the acid-soluble fraction of litter
    REAL(wp)     :: frac_litter_2_water       ! Fraction of incoming material which enters the water-soluble fraction of litter
    REAL(wp)     :: frac_litter_2_ethanol     ! Fraction of incoming material which enters the ethanol-soluble fraction of litter
    REAL(wp)     :: frac_litter_2_nonsoluble  ! Fraction of incoming material which enters the nonsoluble fraction of litter
    REAL(wp)     :: frac_litter_2_humus       ! Fraction of incoming material which enters the humus pool

    ! Yasso parameters (Tuomi '11)
    ! reference decompositation rates for T = 0C and unlimited water availability  in [1/yr]... 
    REAL(wp), PARAMETER   :: ref_decomp_rate_acid       = -0.72_wp    !  ... acid-soluble fraction 
    REAL(wp), PARAMETER   :: ref_decomp_rate_water      = -5.9_wp     !  ... water-soluble fraction 
    REAL(wp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28_wp    !  ... ethanol-soluble fraction 
    REAL(wp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031_wp   !  ... nonsoluble-soluble fraction 
    REAL(wp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016_wp  !  ... humus 
    ! parameters for temperature dependence
    REAL(wp), PARAMETER   :: temp_p1                    = 0.095_wp    ! first  order temperature dependence [1/C]
    REAL(wp), PARAMETER   :: temp_p2                    = -0.0014_wp  ! second order temperature dependence [1/C]
    ! parameter for precipitation dependence
    REAL(wp), PARAMETER   :: precip_p1                  = -1.21_wp    ! first  order precipitation dependence [a/m]
    ! parameters for size dependence 
    REAL(wp), PARAMETER   :: size_p1                    = -1.71_wp    ! first  order size dependence [1/cm]
    REAL(wp), PARAMETER   :: size_p2                    = 0.86_wp     ! second order size dependence [1/cm]
    REAL(wp), PARAMETER   :: size_p3                    = -0.306_wp   ! size dependence power [ ]

    ! Fractions to spit the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere 
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools, 
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux 
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(wp), PARAMETER   :: A_2_W                      =  0.99_wp
    REAL(wp), PARAMETER   :: A_2_E                      =  0.00_wp
    REAL(wp), PARAMETER   :: A_2_N                      =  0.00_wp
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(wp), PARAMETER   :: W_2_A                      = 0.48_wp
    REAL(wp), PARAMETER   :: W_2_E                      = 0.00_wp
    REAL(wp), PARAMETER   :: W_2_N                      = 0.015_wp
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(wp), PARAMETER   :: E_2_A                      = 0.01_wp
    REAL(wp), PARAMETER   :: E_2_W                      = 0.00_wp
    REAL(wp), PARAMETER   :: E_2_N                      = 0.95_wp
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(wp), PARAMETER   :: N_2_A                      = 0.83_wp
    REAL(wp), PARAMETER   :: N_2_W                      = 0.01_wp
    REAL(wp), PARAMETER   :: N_2_E                      = 0.02_wp
    ! fraction of decomposition fluxes of AWEN pools which enters the humus pool [ ]
    REAL(wp), PARAMETER   :: AWEN_2_H                   = 0.0045_wp

    ! climatic drivers
    REAL(wp)     :: precip                  ! 15 runnning mean of precipitation      [m/a]
    REAL(wp)     :: temp                    ! 15 runnning mean of 2m air temperature [C]

    REAL(wp)     :: d_temp                  ! term which accounts for the temperature influence on litter decomposition
    REAL(wp)     :: d_precip                ! term which accounts for the precipitation influence on litter decomposition
    REAL(wp)     :: d_size                  ! term which accounts for the litter size influence on litter decomposition
     
    REAL(wp)     :: redFactor  ! dummy for redFact_Nlimit

    ! time stepping; the parameterization of yasso is done on a annual time
    ! step, thus the fluxes have to be scaled up to annual rates. 
    REAL(wp)     :: dt
    dt = 1.0_wp/days_per_year   

    IF(PRESENT(redFact_Nlimit)) THEN ! nitrogen cycle active & we know already the Nlimitation factor
        redFactor  = redFact_Nlimit
    ELSE
        redFactor  = 1.0_wp
    ENDIF

    ! 0. Preparations
     ! initialize parameters
     frac_litter_2_acid       = Lit_coefV(1)
     frac_litter_2_water      = Lit_coefV(2)
     frac_litter_2_ethanol    = Lit_coefV(3)
     frac_litter_2_nonsoluble = Lit_coefV(4)
     frac_litter_2_humus      = Lit_coefV(5)

     ! initialize yasso pools (AWEN + humus)
     YCpool_acid_ag        = Yasso_io_pools(1)
     YCpool_water_ag       = Yasso_io_pools(2)
     YCpool_ethanol_ag     = Yasso_io_pools(3)
     YCpool_nonsoluble_ag  = Yasso_io_pools(4)

     YCpool_acid_bg        = Yasso_io_pools(5)
     YCpool_water_bg       = Yasso_io_pools(6)
     YCpool_ethanol_bg     = Yasso_io_pools(7)
     YCpool_nonsoluble_bg  = Yasso_io_pools(8)

     YCpool_humus          = Yasso_io_pools(9)

    ! get the sum of the AWEN pools before decomposition
    Pseudo_litter_green = YCpool_acid_ag          &
                        + YCpool_water_ag         &
                        + YCpool_ethanol_ag       &
                        + YCpool_nonsoluble_ag    &
                        + YCpool_acid_bg          &
                        + YCpool_water_bg         &
                        + YCpool_ethanol_bg       &
                        + YCpool_nonsoluble_bg

     ! Change units of climatic forcing variables
     precip = Weather(2)*sec_per_year/1000._wp  ! mm/s -> m/a 
     temp   = Weather(1) - Tmelt                ! K -> C

     ! Calculate annual litter influxes [mol(c)/m2/a]
     Cflx_litter_2_acid       = frac_litter_2_acid        *litter * days_per_year
     Cflx_litter_2_water      = frac_litter_2_water       *litter * days_per_year
     Cflx_litter_2_ethanol    = frac_litter_2_ethanol     *litter * days_per_year
     Cflx_litter_2_nonsoluble = frac_litter_2_nonsoluble  *litter * days_per_year
     Cflx_litter_2_humus      = frac_litter_2_humus       *litter * days_per_year


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 1. Calculate decomposition rates 

     ! Temperature dependence of decomposition
     d_temp   = EXP(temp_p1*temp + temp_p2*temp**2.0_wp)  
     ! Precipitation dependence of decomposition
     d_precip = 1.0_wp - EXP(precip_p1*precip)
     ! Litter size dependence of decomposition -- no effect if WoodlitterSize = 0.0 
     d_size   = MIN(1.0_wp,(1.0_wp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_wp)**size_p3)

     ! decomposition rates accounting for temperature, precipitation, litter size, and nutrient limitation 
     d_acid       =  redFactor * ref_decomp_rate_acid       *d_temp *d_precip *d_size
     d_water      =  redFactor * ref_decomp_rate_water      *d_temp *d_precip *d_size
     d_ethanol    =  redFactor * ref_decomp_rate_ethanol    *d_temp *d_precip *d_size
     d_nonsoluble =  redFactor * ref_decomp_rate_nonsoluble *d_temp *d_precip *d_size

     ! the decomposition of humus is not limited by nutrients; thus no redFactor
     d_humus      =  ref_decomp_rate_humus      *d_temp *d_precip           ! no size effect on humus 

    ! 2. Calculate fluxes 

     ! loss fluxes (negative values): 
     Cflx_from_acid_ag         = YCpool_acid_ag       * d_acid
     Cflx_from_water_ag        = YCpool_water_ag      * d_water
     Cflx_from_ethanol_ag      = YCpool_ethanol_ag    * d_ethanol
     Cflx_from_nonsoluble_ag   = YCpool_nonsoluble_ag * d_nonsoluble

     Cflx_from_acid_bg         = YCpool_acid_bg       * d_acid
     Cflx_from_water_bg        = YCpool_water_bg      * d_water
     Cflx_from_ethanol_bg      = YCpool_ethanol_bg    * d_ethanol
     Cflx_from_nonsoluble_bg   = YCpool_nonsoluble_bg * d_nonsoluble

     Cflx_from_humus           = YCpool_humus         * d_humus

     ! gain fluxes (positive values): 
     ! fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why) 
     Cflx_2_acid_ag           = ABS(Cflx_from_water_ag      * W_2_A  &    ! returns positive fluxes
                                  + Cflx_from_ethanol_ag    * E_2_A  &
                                  + Cflx_from_nonsoluble_ag * N_2_A) 

     Cflx_2_water_ag          = ABS(Cflx_from_acid_ag       * A_2_W  &
                                  + Cflx_from_ethanol_ag    * E_2_W  &
                                  + Cflx_from_nonsoluble_ag * N_2_W) 

     Cflx_2_ethanol_ag        = ABS(Cflx_from_acid_ag       * A_2_E  &
                                  + Cflx_from_water_ag      * W_2_E  &
                                  + Cflx_from_nonsoluble_ag * N_2_E) 

     Cflx_2_nonsoluble_ag     = ABS(Cflx_from_acid_ag       * A_2_N  &
                                  + Cflx_from_water_ag      * W_2_N  &
                                  + Cflx_from_ethanol_ag    * E_2_N) 

     Cflx_2_acid_bg           = ABS(Cflx_from_water_bg      * W_2_A  &    
                                  + Cflx_from_ethanol_bg    * E_2_A  &
                                  + Cflx_from_nonsoluble_bg * N_2_A) 

     Cflx_2_water_bg          = ABS(Cflx_from_acid_bg       * A_2_W  &
                                  + Cflx_from_ethanol_bg    * E_2_W  &
                                  + Cflx_from_nonsoluble_bg * N_2_W) 

     Cflx_2_ethanol_bg        = ABS(Cflx_from_acid_bg       * A_2_E  &
                                  + Cflx_from_water_bg      * W_2_E  &
                                  + Cflx_from_nonsoluble_bg * N_2_E) 

     Cflx_2_nonsoluble_bg     = ABS(Cflx_from_acid_bg       * A_2_N  &
                                  + Cflx_from_water_bg      * W_2_N  &
                                  + Cflx_from_ethanol_bg    * E_2_N) 

     Cflx_2_humus          = ABS(Cflx_from_acid_ag                &
                               + Cflx_from_water_ag               &
                               + Cflx_from_ethanol_ag             &
                               + Cflx_from_nonsoluble_ag          &
                               + Cflx_from_acid_bg                &
                               + Cflx_from_water_bg               &
                               + Cflx_from_ethanol_bg             &
                               + Cflx_from_nonsoluble_bg          &
                               ) * AWEN_2_H  

     Cflx_2_humusAG        = ABS(Cflx_from_acid_ag                &
                               + Cflx_from_water_ag               &
                               + Cflx_from_ethanol_ag             &
                               + Cflx_from_nonsoluble_ag          &
                               ) * AWEN_2_H

     Cflx_2_humusBG        = ABS(Cflx_from_acid_bg                &
                               + Cflx_from_water_bg               &
                               + Cflx_from_ethanol_bg             &
                               + Cflx_from_nonsoluble_bg          &
                               ) * AWEN_2_H

     ! the remaining fractions of the loss fluxes enter the atmosphere as respiration
     soilResp_rateYasso    =  (Cflx_from_acid_ag + Cflx_from_acid_bg)               &
                              * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                              + (Cflx_from_water_ag + Cflx_from_water_bg)           &
                              * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                              + (Cflx_from_ethanol_ag + Cflx_from_ethanol_bg)       &
                              * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                              + (Cflx_from_nonsoluble_ag + Cflx_from_nonsoluble_bg) &
                              * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)         &
                              + Cflx_from_humus 

     ! litter ag & bg respiration (needed for N cycle)                       
     soilResp_rateLitterAG   =  (Cflx_from_acid_ag )                                &
                              * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                              + (Cflx_from_water_ag )                               &
                              * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                              + (Cflx_from_ethanol_ag )                             &
                              * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                              + (Cflx_from_nonsoluble_ag )                          &
                              * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)         

     soilResp_rateLitterBG   =  (Cflx_from_acid_bg )                                &
                              * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                              + (Cflx_from_water_bg )                               &
                              * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                              + (Cflx_from_ethanol_bg )                             &
                              * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                              + (Cflx_from_nonsoluble_bg )                          &
                              * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)

   ! 3. update Yasso pools
     YCpool_acid_ag           = MAX(0.0_wp,                                         &
                              YCpool_acid_ag                                        & ! old pool   
                              + (Cflx_from_acid_ag                                  & ! incoming flux from AWEN pools
                                + Cflx_2_acid_ag                                    & ! outgoing flux
                                + frac_aboveground * Cflx_litter_2_acid) *dt)         ! incoming flux from litter
     YCpool_water_ag          = MAX(0.0_wp,                                         & ! and so on .....
                              YCpool_water_ag                                       &
                              + (Cflx_from_water_ag                                 &
                                + Cflx_2_water_ag                                   &
                                + frac_aboveground * Cflx_litter_2_water) *dt)
     YCpool_ethanol_ag        = MAX(0.0_wp,                                         &
                              YCpool_ethanol_ag                                     &
                              + (Cflx_from_ethanol_ag                               &
                                + Cflx_2_ethanol_ag                                 &
                                + frac_aboveground * Cflx_litter_2_ethanol) *dt)
     YCpool_nonsoluble_ag     = MAX(0.0_wp,                                         &
                              YCpool_nonsoluble_ag                                  &
                              + (Cflx_from_nonsoluble_ag                            &
                                + Cflx_2_nonsoluble_ag                              &
                                + frac_aboveground * Cflx_litter_2_nonsoluble) *dt)
    

     YCpool_acid_bg           = MAX(0.0_wp,                                         &
                              YCpool_acid_bg                                        & ! old pool   
                              + (Cflx_from_acid_bg                                  & ! incoming flux from AWEN pools
                                + Cflx_2_acid_bg                                    & ! outgoing flux
                                + (1._wp - frac_aboveground) *                      & ! 
                                  Cflx_litter_2_acid) *dt)                            ! incoming flux from litter
     YCpool_water_bg          = MAX(0.0_wp,                                         & ! ...
                              YCpool_water_bg                                       &
                              + (Cflx_from_water_bg                                 &
                                + Cflx_2_water_bg                                   &
                                + (1._wp - frac_aboveground) *                      &
                                  Cflx_litter_2_water) *dt                          &
                              + NPP_2_rootExudates)                                    ! exudates are carbohydrates only and belowground
     YCpool_ethanol_bg        = MAX(0.0_wp,                                         &
                              YCpool_ethanol_bg                                     &
                              + (Cflx_from_ethanol_bg                               &
                                + Cflx_2_ethanol_bg                                 &
                                + (1._wp - frac_aboveground) *                      &
                                     Cflx_litter_2_ethanol) *dt)
     YCpool_nonsoluble_bg     = MAX(0.0_wp,                                         &
                              YCpool_nonsoluble_bg                                  &
                              + (Cflx_from_nonsoluble_bg                            &
                                + Cflx_2_nonsoluble_bg                              &
                                + (1._wp - frac_aboveground) *                      &
                                  Cflx_litter_2_nonsoluble) *dt)

     YCpool_humus          = MAX(0.0_wp,                            &
                              YCpool_humus                          &
                              + (Cflx_from_humus                    &  
                                + Cflx_2_humus                      &
                                + Cflx_litter_2_humus) *dt)

     ! compute d_litter_green; this is given by the change in AWEN pools
     !          Flx                 Pool(t)                   Pool(t+1)
     !    d_P = ----     &   Flx = ----------   --->   d_P = -----------  - 1
     !          Pool(t)             Pool(t+1)                 Pool(t)

     IF (Pseudo_litter_green .GT. 0.0_wp) THEN
        d_litter_green = (YCpool_acid_ag            &
                          + YCpool_water_ag         &
                          + YCpool_ethanol_ag       &
                          + YCpool_nonsoluble_ag    &
                          + YCpool_acid_bg          &
                          + YCpool_water_bg         &
                          + YCpool_ethanol_bg       &
                          + YCpool_nonsoluble_bg)   &
                          / Pseudo_litter_green - 1._wp
      ELSE
         d_litter_green = 0.0_wp
      ENDIF

   ! 4.1 Write pools into the output array
    Yasso_out(1)           = YCpool_acid_ag 
    Yasso_out(2)           = YCpool_water_ag
    Yasso_out(3)           = YCpool_ethanol_ag
    Yasso_out(4)           = YCpool_nonsoluble_ag
    Yasso_out(5)           = YCpool_acid_bg 
    Yasso_out(6)           = YCpool_water_bg
    Yasso_out(7)           = YCpool_ethanol_bg
    Yasso_out(8)           = YCpool_nonsoluble_bg
    Yasso_out(9)           = YCPool_humus
   ! 4.2 Write respiration into the output array
    Yasso_out(10)           = soilResp_rateYasso * dt ! convert back to daily

    IF (.NOT.PRESENT(redFact_Nlimit)) THEN
       ! the routine might have been called to diagnose the nutrient demand, so we give out 
       ! 1. flux for immobilisation: this is the sum of all carbon coming from the non-humus pools to humus
       Yasso_out(11)           = Cflx_2_humusAG * dt ! daily ! the flux is positive
       Yasso_out(12)           = Cflx_2_humusBG * dt ! daily ! the flux is positive
       ! 2. flux for mineralisation:  this is the flux from humus to atmosphere (Cflx_slow_2_atmos)
       Yasso_out(13)           =  - Cflx_from_humus * dt ! daily ! sign CHECKED!
       ! ... as well as all fluxes from AWEN pools to atmosphere
       Yasso_out(14)           = -soilResp_rateLitterAG * dt !daily ! the flux is negative
       Yasso_out(15)           = -soilResp_rateLitterBG * dt !daily ! the flux is negative
       ! 3. decomposition rate for Npool litter green 
       Yasso_out(16)           = d_litter_green
    ENDIF

  END SUBROUTINE yasso


  ! --- update_delta_carbon() ---------------------------------------------------------------------------------------------------
  !
  ! updates delta carbon on area shift for a single pool 
  !
  elemental subroutine update_delta_carbon( &
    & org_Cpool,                  & ! in
    & delta_Cpool_source,         & ! inout
    & delta_Cpool_sink,           & ! inout
    & delta_fract,                & ! in
    & prev_fract)                   ! in

    real(wp),intent(in)    ::  org_Cpool, delta_fract, prev_fract
    real(wp),intent(inout) ::  delta_Cpool_source, delta_Cpool_sink

    !TODO-JN: probably needs to be adapted with absolute areas, then also the prev_fract 
    !         will be required + the test for zero 
!    IF (delta_fract .GT. 0) THEN
        delta_Cpool_source = delta_Cpool_source - (org_Cpool * delta_fract)
        delta_Cpool_sink = delta_Cpool_sink + (org_Cpool * delta_fract)
!    ENDIF

  end subroutine update_delta_carbon


  ! --- update_state_reset_delta() ---------------------------------------------------------------------------------------------------
  !
  ! add delta carbon to carbon state and reset to zero 
  !
  elemental subroutine update_state_reset_delta( &
    & org_Cpool,                  & ! inout
    & delta_Cpool,                & ! inout
    & delta_fract,                & ! in
    & prev_fract)                   ! in

    real(wp),intent(in)    ::  delta_fract, prev_fract
    real(wp),intent(inout) ::  org_Cpool, delta_Cpool

    !TODO-JN: probably needs to be adapted with absolute areas

    ! TODO-JN: do they still need smallest possible fraction 1e-10? 
    ! (Then the following should not happen and other Assertions and MAX statements will be required!)
    ! IF prev_fract + delta_fract .EQ. 0._wp org_Cpool is not touched
    ! i.e. the densities are not altered -- thus empty tiles can have non zero densities
    ! -- if smallest possible fraction    
    IF( prev_fract + delta_fract .GT. 0._wp) THEN
      ! i.e. area gained
      org_Cpool = ((org_Cpool * prev_fract) + delta_Cpool) / (prev_fract + delta_fract)
    ENDIF
    delta_Cpool = 0._wp

  end subroutine update_state_reset_delta


  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2 
  !
  elemental subroutine get_per_tile( &
    & Cpool_per_tile,                & ! in
    & Cpool_per_canopy,              & ! in
    & veg_frac_correction,           & ! in
    & fract_fpc_max)                   ! in

    real(wp),intent(in) ::  fract_fpc_max, veg_frac_correction
    real(wp),intent(in) ::  Cpool_per_canopy
    real(wp),intent(out)::  Cpool_per_tile

    Cpool_per_tile = Cpool_per_canopy * veg_frac_correction * fract_fpc_max

  end subroutine get_per_tile

#endif
END MODULE mo_carbon_process
