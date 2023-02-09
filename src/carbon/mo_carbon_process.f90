!>
!! @brief Contains the routines for the carbon processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!! Julia Nabel (2018-01-10): l_forestRegrowth flag -> allow for changing max_LAI
!! Julia Nabel (2021-05-28): start of some refactoring (more advisable)

!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_carbon_process
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround, &
                               & days_per_year, sec_per_year, sec_per_day,        &
                               & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, &
                               & i_lctlib_nonsoluble, i_lctlib_humus
  USE mo_carbon_process_shared, ONLY: yasso

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_Cpools, relocate_carbon_fire, relocate_carbon_damage

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_process'

CONTAINS

  ! Calculations to allocate NPP to the plant C pools.
  ELEMENTAL SUBROUTINE calc_Cpools( &
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
    & soilResp_rate,              & ! out
    & NPP_flux_correction,        & ! out
    & excess_NPP,                 & ! out
    & root_exudates,              & ! out
    & cflux_c_greenwood_2_litter, & ! out
    & cflux_herbivory,            & ! out
    & cflux_herb_2_littergreen,   & ! out
    & cflux_herb_2_atm,           & ! out
    & NPP_act,                    & ! out
    !
    ! & fract_litter_wood_new,      & ! inout, optional; R: only for spitfire
    ! Yasso variables
    & temp2_30d,                  & ! in
    & precip_30d,                 & ! in
    !
    & c_acid_ag1,                 & ! inout
    & c_water_ag1,                & ! inout
    & c_ethanol_ag1,              & ! inout
    & c_nonsoluble_ag1,           & ! inout
    & c_acid_bg1,                 & ! inout
    & c_water_bg1,                & ! inout
    & c_ethanol_bg1,              & ! inout
    & c_nonsoluble_bg1,           & ! inout
    & c_humus_1,                  & ! inout
    & c_acid_ag2,                 & ! inout
    & c_water_ag2,                & ! inout
    & c_ethanol_ag2,              & ! inout
    & c_nonsoluble_ag2,           & ! inout
    & c_acid_bg2,                 & ! inout
    & c_water_bg2,                & ! inout
    & c_ethanol_bg2,              & ! inout
    & c_nonsoluble_bg2,           & ! inout
    & c_humus_2,                  & ! inout
    !
    & LeafLit_coef_acid,          & ! in
    & LeafLit_coef_water,         & ! in
    & LeafLit_coef_ethanol,       & ! in
    & LeafLit_coef_nonsoluble,    & ! in
    & LeafLit_coef_humus,         & ! in
    & WoodLit_coef_acid,          & ! in
    & WoodLit_coef_water,         & ! in
    & WoodLit_coef_ethanol,       & ! in
    & WoodLit_coef_nonsoluble,    & ! in
    & WoodLit_coef_humus,         & ! in
    & WoodLitterSize,             & ! in
    !
    & c_decomp_humus_1_sum,       & ! out, optional, necessary with yasso
    & c_decomp_humus_2_sum,       & ! out, optional, necessary with yasso
    & c_into_humus_1_sum,         & ! out, optional, necessary with yasso
    & c_into_humus_2_sum          & ! out, optional, necessary with yasso
    & )

    USE mo_jsb_physical_constants, ONLY : tmelt     ! tmelt =273.15_wp

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
    !
    REAL(wp), INTENT(out)   :: soilResp_rate          ! Soil (=heterotrophic) respiration rate  [mol(C)/m^2 s] Note: this is a loss,
                                                     !    therefore negative!
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
    !
    ! added for thonicke fire algorithm (spitfire) needed to correctly split the carbon pools into the different fuel classes
    !real(wp),optional,intent(inout) :: fract_litter_wood_new   ! new fraction in above ground wood litter pool

    ! Meteorology for Yasso
    REAL(wp), INTENT(in)    :: temp2_30d              ! 30 day mean temperature
    REAL(wp), INTENT(in)    :: precip_30d             ! 30 day mean precipitation

    ! Yasso pools
    ! Size class 1; green litter
    !   Aboveground
    REAL(wp), INTENT(inout) :: c_acid_ag1        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag1       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag1     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag1  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground
    REAL(wp), INTENT(inout) :: c_acid_bg1        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg1       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg1     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_1         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2; woody litter
    !   Aboveground
    REAL(wp), INTENT(inout) :: c_acid_ag2        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag2       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag2     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag2  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground
    REAL(wp), INTENT(inout) :: c_acid_bg2        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg2       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg2     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_2         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Parameters for Yasso
       ! Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(wp), INTENT(in)    :: LeafLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_water         ! fraction going to the water soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_humus         ! fraction going to the humus pool
       ! Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(wp), INTENT(in)    :: WoodLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_water         ! fraction going to the water soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_humus         ! fraction going to the humus pool
    REAL(wp), INTENT(in)    :: WoodLitterSize             ! size of coarse debris

    !
    REAL(wp), INTENT(inout) :: c_decomp_humus_1_sum   ! annual cflux sum of humus decomposed (leaf) [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_decomp_humus_2_sum   ! annual cflux sum of humus decomposed (wood) [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_into_humus_1_sum     ! annual cflux sum into humus (leaf) [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_into_humus_2_sum     ! annual cflux sum into humus (wood) [mol(C)/m^2(canopy)]

    ! locals
    LOGICAL  :: is_crop                 ! logical crop mask

    REAL(wp) :: c_green_max,c_reserve_optimal
    REAL(wp) :: NPP_2_greenPool, NPP_2_woodPool, NPP_2_reservePool, NPP_2_rootExudates, Green_2_herbivory
    REAL(wp) :: C_2_litter_greenPool, C_2_litter_woodPool
    REAL(wp) :: excess_carbon
    REAL(wp) :: leaf_shedding

    REAL(wp) :: c_green_pot, c_woods_pot         !! values of potential carbon pools i.e. before accounting for N-limitation
    REAL(wp) :: c_reserve_pot

    REAL(wp) :: Cflx_NPP_2_green_pot,Cflx_NPP_2_wood_pot,Cflx_NPP_2_reserve_pot
    REAL(wp) :: Cflx_2_crop_harvest, Cflx_crop_harvest_2_atmos
    REAL(wp) :: Cflx_faeces_2_LG, Cflx_faeces_2_atm
    ! R: only for spitfire
    !REAL(wp) :: litter_wood_new, litter_wood_total

    ! local variables for Yasso
    REAL(wp), DIMENSION(18) :: Yasso_io_pools
    REAL(wp), DIMENSION(2)  :: Weather
    REAL(wp), DIMENSION(18) :: Yasso_out
    REAL(wp), DIMENSION(5)  :: LeafLit_coefV
    REAL(wp), DIMENSION(5)  :: WoodLit_coefV
    REAL(wp)                :: leafLitter      ! amount of leafLitter
    REAL(wp)                :: sizeZero

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

    !LOGICAL  :: with_Nitrogen        !! .FALSE.: only Carbon pools are updated; .TRUE.: in addition also Nitrogen pools are updated
    !LOGICAL  :: with_yasso           !! .FALSE.: cbalance model for litter and soil carbon is used; .TRUE.: yasso model is used
    !REAL(wp), PARAMETER :: N2O_ef_grazing = 0._wp  !! N2O grazing switched off; reasonable value for the emission factor: 0.0125_wp


    ! R: Normally parameters should be in mo_carbon_constants, but the following parameters are used only in this subroutine!
    !    Therefore I put them here. Is there any not structural advantage to put them into mo_carbon_constants?
    ! Parameter for Subroutine calc_Cpools
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
    REAL(wp), PARAMETER :: N2O_rate_nitrification = 0.001_wp        ! JSBACH assumes lower bound of estimation in the literature as:
                                                                    !   (0.1%) due to nitrification
    REAL(wp), PARAMETER :: N2O_rate_denitrification = 0.00125_wp    ! 0.0025_wp  ! & (0.2%) due to denitrification
                                                                    ! actually 0.2% of denitrified N, not of N inputs
                                                                    ! assumption: 50% loss of NO3 by leaching
    REAL(wp), PARAMETER :: N2O_rate               = 0.005_wp        ! JSBACH assumes 0.5-1%
                                                                    !    N2O emissions rate per day (<0.1%-0.2% nitrification );
                                                                    !    0.2-4.7% due to denit in Xu Ri, 2008 * therein
    REAL(wp), PARAMETER :: sminn_NH4_fraction     = 0.4_wp          ! soil mineral NH4 is 40% of SMINN pool and rest in the form of
                                                                    !  NO3. Ref: Xu Ri, 2008

    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    soilResp_rate            = 0.0_wp
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

    ! with yasso:
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

    Weather(1) = temp2_30d
    Weather(2) = precip_30d

    LeafLit_coefV(i_lctlib_acid)       = LeafLit_coef_acid
    LeafLit_coefV(i_lctlib_water)      = LeafLit_coef_water
    LeafLit_coefV(i_lctlib_ethanol)    = LeafLit_coef_ethanol
    LeafLit_coefV(i_lctlib_nonsoluble) = LeafLit_coef_nonsoluble
    LeafLit_coefV(i_lctlib_humus)      = LeafLit_coef_humus
    WoodLit_coefV(i_lctlib_acid)       = WoodLit_coef_acid
    WoodLit_coefV(i_lctlib_water)      = WoodLit_coef_water
    WoodLit_coefV(i_lctlib_ethanol)    = WoodLit_coef_ethanol
    WoodLit_coefV(i_lctlib_nonsoluble) = WoodLit_coef_nonsoluble
    WoodLit_coefV(i_lctlib_humus)      = WoodLit_coef_humus

    Yasso_out(1:18) = 0.0_wp           ! array to store the Yasso pools/fluxes
    leafLitter      = 0.0_wp           ! amount of leafLitter
    sizeZero        = 0.0_wp

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
    excess_carbon = max(0.0_wp,c_green - c_green_max)     !! If by the reduction of LAI the maximum value of the green
                                                                  !!    pool is still smaller than the current value, there is
    C_2_litter_greenPool = C_2_litter_greenPool + excess_carbon   !!    excess C that has also to be shedded to the green litter
    c_green = c_green - excess_carbon                     !!    pool and also subtracted from the green pool

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
      c_woods_pot = c_woods_pot - excess_carbon                   !! .. so that the excess carbon has once more to be
                                                                  !! .. subtracted
      NPP_2_greenPool = NPP_2_greenPool + excess_carbon           !! .. and instead made available to the green pool
      NPP_2_woodPool = NPP_2_woodPool - excess_carbon             !! .. and the actual amount of NPP put to the wood pool
                                                                  !! .. is less
      !! Growth of Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
      IF (c_reserve_pot < c_reserve_optimal) THEN                 !! .. If the reserve pool is smaller than optimal,
        c_reserve_pot = c_reserve_pot +  NPP_2_reservePool        !! .... it is filled by the available NPP
        excess_carbon =                                         & !! .... Thereby it may happen that it gets larger than
             max(0.0_wp,c_reserve_pot - c_reserve_optimal)        !! .... optimal so that there is excess carbon
        c_reserve_pot = c_reserve_pot - excess_carbon             !! .... that needs not be taken up
        NPP_2_greenPool = NPP_2_greenPool + excess_carbon         !! .... but can better be used to increase the green
        NPP_2_reservePool = NPP_2_reservePool - excess_carbon     !! .... pool and the actual amount of NPP put to the
                                                                  !! .... reserve pool is less.
      ELSE                                                        !! .... Otherwise (reserve pool is larger than optimal)
        NPP_2_greenPool = NPP_2_greenPool + NPP_2_reservePool     !! .... all NPP is left for the green pool
        NPP_2_reservePool = 0.0_wp                                !! .... so that nothing is stored in the reserve pool.
      END IF

      !! Growth of Green Pool (leaves and fine roots): (in case of too much NPP, try to put it into the reserve pool).

      c_green_pot = c_green_pot + NPP_2_greenPool                 !! Green pool is filled by the available NPP.
      excess_carbon = max(0.0_wp,c_green_pot - c_green_max)       !! .. Thereby it may get larger than appropriate for
                                                                  !!    current LAI.
      NPP_2_greenPool = NPP_2_greenPool - excess_carbon           !! .. Hence the actual amount of NPP put to the green
                                                                  !!    pool is less.
      c_green_pot = c_green_pot - excess_carbon                   !! .. This excess carbon needs not be taken up but
      IF (c_reserve_pot < c_reserve_optimal) THEN                 !! .... if the reserve pool is smaller than optimal
        c_reserve_pot = c_reserve_pot + excess_carbon             !! .... it is tried to put the carbon there,
        NPP_2_reservePool = NPP_2_reservePool + excess_carbon     !! .... which means that additional NPP is put to the
                                                                  !!      reserve pool.
        excess_carbon =                                        &  !! .... Thereby it may happen that the reserve pool
             max(0.0_wp,c_reserve_pot - c_reserve_optimal)        !!      increases beyond the optimal value.
        c_reserve_pot = c_reserve_pot - excess_carbon             !! .... In that case the excess carbon is once more
                                                                  !!      removed from the reserve pool
        NPP_2_reservePool = NPP_2_reservePool - excess_carbon     !! .... so that the actual amount of NPP put into the
                                                                  !!      reserve pool is less.
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
    c_reserve     = c_reserve     +  Cflx_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N

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
    soilResp_rate = 0._wp

    !! call yasso for leaf litter; this time with nutrient effects
    CALL yasso(Yasso_io_pools(1:9), Weather, leafLitter, LeafLit_coefV,               &
      &        sizeZero, Yasso_out, fract_green_aboveGround, NPP_2_rootExudates,redFact_Nlimitation)

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
    soilResp_rate = Yasso_out(10) / sec_per_day

    !! call yasso for woody litter; this time with nutrient effects
    CALL yasso(Yasso_io_pools(10:18), Weather, C_2_litter_woodPool, WoodLit_coefV,       &
      & WoodLitterSize, Yasso_out, fract_wood_aboveGround, 0.0_wp, redFact_Nlimitation) ! no exudates

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
    soilResp_rate = soilResp_rate  + Yasso_out(10) / sec_per_day

    c_decomp_humus_2_sum = c_decomp_humus_2_sum + Yasso_out(17)
    c_into_humus_2_sum   = c_into_humus_2_sum + Yasso_out(18)


    ! always:
    !! ============== Decay of crop harvest =============

    Cflx_crop_harvest_2_atmos = c_crop_harvest / tau_c_crop_harvest
    c_crop_harvest = c_crop_harvest - Cflx_crop_harvest_2_atmos
    soilResp_rate = soilResp_rate - Cflx_crop_harvest_2_atmos / sec_per_day
    !soilResp_rate_pot = soilResp_rate_pot - Cflx_crop_harvest_2_atmos / sec_per_day

    !! ============== COMPUTE Herbivory RESPIRATION =====

    cflux_herb_2_atm = -(Cflx_faeces_2_atm) /sec_per_day

    !! ============== Litter and herbivory diagnostics (flux from day-1 to s-1)
    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter /sec_per_day
    cflux_herbivory    = Green_2_herbivory / sec_per_day
    cflux_herb_2_littergreen = Cflx_faeces_2_LG / sec_per_day   ! DSG: goes into yasso



    ! R: only for spitfire:
    ! fraction of new wood litter, assume old and new litter are respired in the same way, respiration does not change the fraction
    !        if (present(fract_litter_wood_new)) then
    !           if (with_yasso) then
    ! litter_wood_total = c_water_ag2 + c_acid_ag2 + c_ethanol_ag2 + c_nonsoluble_ag2
    !           else
    !             litter_wood_total = c_litter_wood_ag
    !           endif
    !           if ((litter_wood_total >  EPSILON(1._wp)) .AND. (litter_wood_new < litter_wood_total)) then
    !              fract_litter_wood_new = litter_wood_new / litter_wood_total
    !           else
    !              fract_litter_wood_new = 1._wp    ! fract_litter_wood_new=1 means fuel fractions will be reset to original values
    !           endif
    !        endif

  END SUBROUTINE calc_Cpools




  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of fire for the carbon pools. More precisely:
  ! It is assumed that for the burned area the carbon from the above ground litter pools (c_litter_green_ag,
  ! c_litter_wood_ag) is completely released to the atmosphere and the carbon from the living plant pools
  ! (c_green_st, c_reserve_st, c_woods_st) is partly released to the atmosphere and partly relocated
  ! into the litter pools.
  !
  SUBROUTINE relocate_carbon_fire( &
    & burned_fract,                      & ! in from disturbed_frac
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
    REAL(wp), INTENT(in)    :: burned_fract(:)              ! Fraction of the vegetated area of each tile burned till the last
                                                            ! call of this routine
    REAL(wp), INTENT(in)    :: fire_fract_wood_2_atmos      ! Fraction of above ground wood immediately emitted to atmosphere
                                                            ! by fire
    REAL(wp), INTENT(in)    :: LeafLit_coef(:)              ! Factor to spread non woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(in)    :: WoodLit_coef(:)              ! Factor to spread woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(inout) :: c_green(:)                   ! Value of green carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_woods(:)                   ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_reserve(:)                 ! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(out)   :: co2flux_fire_all_2_atm(:)    ! CO2 immediately released by fire per tile [kg(CO2)m-2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_green_2_soil(:)   ! Amount of carbon relocated by wind and fire damage
                                                            ! .. to the green litter pools [mol(C)/m^2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   ! Amount of carbon relocated by wind damage and fire
                                                            ! .. to the wood litter pools [mol(C)/m^2(canopy)s-1]
    !   YASSO
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    REAL(wp), INTENT(inout) :: c_acid_ag1(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag1(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag1(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag1(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(wp), INTENT(inout) :: c_acid_bg1(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg1(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg1(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_1(:)           ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    REAL(wp), INTENT(inout) :: c_acid_ag2(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag2(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag2(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag2(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(wp), INTENT(inout) :: c_acid_bg2(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg2(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg2(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_2(:)           ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp), DIMENSION(SIZE(burned_fract))  :: non_woody_litter, fraction ! Helper
    INTEGER                                  :: ind, nc

    nc = SIZE(burned_fract)
    non_woody_litter(:) = c_green(:) + c_reserve(:)

    !! diagnose amount of carbon released from the green and reserve pool to the green litter pools
    !! diagnose amount of carbon released from the wood pool to the woody litter pools
    cflux_dist_green_2_soil(:) = cflux_dist_green_2_soil(:) + &
      & (c_green(:) + c_reserve(:)) * (1._wp - fract_green_aboveGround) * burned_fract(:) / sec_per_day
    cflux_dist_woods_2_soil(:)  = cflux_dist_woods_2_soil(:) + &
      & c_woods(:) * (1._wp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * burned_fract(:) / sec_per_day

    !! determine amount of carbon released from the wood pool, living tissue pools, and litter pools to the atmosphere
    co2flux_fire_all_2_atm(:) = ( &
        ! take living parts of plants above ground:
      & non_woody_litter(:) * fract_green_aboveGround                      &
        ! now add the YASSO pools above ground:
      & + c_water_ag1(:) + c_acid_ag1(:) + c_ethanol_ag1(:)                &
      & + c_nonsoluble_ag1(:)                                              &
      & + c_water_ag2(:) + c_acid_ag2(:) + c_ethanol_ag2(:)                &
      & + c_nonsoluble_ag2(:)                                              &
        ! now add the fire affected fraction of wood pools above ground:
      & + c_woods(:) * fire_fract_wood_2_atmos * fract_wood_aboveGround    &
      & ) * burned_fract(:) * molarMassCO2_kg / sec_per_day

    ! Lower down density of aboveground litter pools
    fraction(:) = 1._wp - burned_fract(:)
    c_acid_ag1      (:) = fraction(:) * c_acid_ag1       (:)
    c_water_ag1     (:) = fraction(:) * c_water_ag1      (:)
    c_ethanol_ag1   (:) = fraction(:) * c_ethanol_ag1    (:)
    c_nonsoluble_ag1(:) = fraction(:) * c_nonsoluble_ag1 (:)

    c_acid_ag2      (:) = fraction(:) * c_acid_ag2       (:)
    c_water_ag2     (:) = fraction(:) * c_water_ag2      (:)
    c_ethanol_ag2   (:) = fraction(:) * c_ethanol_ag2    (:)
    c_nonsoluble_ag2(:) = fraction(:) * c_nonsoluble_ag2 (:)

    ! Add green and reserve vegetation carbon to belowground green litter pools
    fraction(:) = (1._wp-fract_green_aboveGround) * burned_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_bg1(ind), c_water_bg1(ind), c_ethanol_bg1(ind), c_nonsoluble_bg1(ind), c_humus_1(ind), &
        &                          non_woody_litter(ind), fraction(ind), LeafLit_coef)
    ENDDO

    ! Add the remaings of burned wood:
    ! AG
    fraction(:) = (1._wp-fire_fract_wood_2_atmos) * fract_wood_aboveGround * burned_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_ag2(ind), c_water_ag2(ind), c_ethanol_ag2(ind), c_nonsoluble_ag2(ind), c_humus_2(ind), &
        &                          c_woods(ind), fraction(ind), WoodLit_coef)
    ENDDO
    ! BG
    fraction(:) = (1._wp-fract_wood_aboveGround) * burned_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_bg2(ind), c_water_bg2(ind), c_ethanol_bg2(ind), c_nonsoluble_bg2(ind), c_humus_2(ind), &
        &                          c_woods(ind), fraction(ind), WoodLit_coef)
    ENDDO

    ! Lower down the carbon density of living plant pools
    fraction(:) = 1._wp - burned_fract(:)
    c_green(:)   = c_green(:)   * fraction(:)
    c_woods(:)   = c_woods(:)   * fraction(:)
    c_reserve(:) = c_reserve(:) * fraction(:)

  END SUBROUTINE relocate_carbon_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (c_green_st, c_reserve_st and
  ! c_woods_st) is partly relocated into the litter pools.
  !
  SUBROUTINE relocate_carbon_damage(  &
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
    real(wp),intent(in)    :: damaged_fract(:)         ! Fraction of the vegetated area of each tile damaged since the last
                                                       ! call of this routine
    real(wp),intent(in)    :: LeafLit_coef(:)          ! Factor to spread non woody litter into yasso pools [ ]
    real(wp),intent(in)    :: WoodLit_coef(:)          ! Factor to spread woody litterfall into yasso pools [ ]
    real(wp),intent(inout) :: c_green(:)               ! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_woods(:)               ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_reserve(:)             ! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(wp),intent(inout) :: cflux_dist_green_2_soil(:)          ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
    real(wp),intent(inout) :: cflux_dist_woods_2_soil(:)          ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
                                                                  ! .. to the green litter pools [mol(C)/m^2(canopy)s-1]

    !   YASSO
    ! Size class 1: green
    !   above ground C pools
    real(wp),intent(inout) :: c_acid_ag1(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_ag1(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_ag1(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_ag1(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: c_acid_bg1(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_bg1(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_bg1(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_bg1(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_humus_1(:)         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2: wood
    !   above ground C pools
    real(wp),intent(inout) :: c_acid_ag2(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_ag2(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_ag2(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_ag2(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: c_acid_bg2(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_bg2(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_bg2(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_bg2(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_humus_2(:)         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp), DIMENSION(SIZE(damaged_fract)) :: non_woody_litter, fraction  ! Helper
    INTEGER                                  :: ind, nc

    ! Initializations
    nc = SIZE(damaged_fract)
    non_woody_litter(:) = c_green(:) + c_reserve(:)

    !! transfer carbon from the wood pool to the wood litter pools nd
    !! transfer carbon from the green and reserve pool to the green litter pool
    cflux_dist_woods_2_soil(:) = cflux_dist_woods_2_soil(:) + c_woods(:) * damaged_fract(:) / sec_per_day
    cflux_dist_green_2_soil(:) = cflux_dist_green_2_soil(:) + non_woody_litter(:) * damaged_fract(:) / sec_per_day

    ! Add damaged green and reserve
    ! AG
    fraction(:) = fract_green_aboveGround * damaged_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_ag1(ind), c_water_ag1(ind), c_ethanol_ag1(ind), c_nonsoluble_ag1(ind), c_humus_1(ind), &
        &                          non_woody_litter(ind), fraction(ind), LeafLit_coef)
    ENDDO
    ! BG
    fraction(:) = (1._wp - fract_green_aboveGround) * damaged_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_bg1(ind), c_water_bg1(ind), c_ethanol_bg1(ind), c_nonsoluble_bg1(ind), c_humus_1(ind), &
        &                          non_woody_litter(ind), fraction(ind), LeafLit_coef)
    ENDDO

    ! Add damaged wood
    ! AG
    fraction(:) = fract_wood_aboveGround  * damaged_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_ag2(ind), c_water_ag2(ind), c_ethanol_ag2(ind), c_nonsoluble_ag2(ind), c_humus_2(ind), &
        &                          c_woods(ind), fraction(ind), WoodLit_coef)
    ENDDO
    ! BG
    fraction(:) = (1._wp - fract_wood_aboveGround)  * damaged_fract(:)
    DO ind = 1,nc
      CALL distribute_yasso_litter(c_acid_bg2(ind), c_water_bg2(ind), c_ethanol_bg2(ind), c_nonsoluble_bg2(ind), c_humus_2(ind), &
        &                          c_woods(ind), fraction(ind), WoodLit_coef)
    ENDDO

    ! Lower down the carbon density of living plant pools
    fraction(:) = 1._wp - damaged_fract(:)
    c_green(:)   = c_green(:)   * fraction(:)
    c_woods(:)   = c_woods(:)   * fraction(:)
    c_reserve(:) = c_reserve(:) * fraction(:)

  END SUBROUTINE relocate_carbon_damage

#endif
END MODULE mo_carbon_process
