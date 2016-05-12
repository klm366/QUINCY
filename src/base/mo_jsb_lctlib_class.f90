! NOTE: This module will be replaced!!!
!>
!! @brief Contains the memory class for landcover type library (lctlib).
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-08)
!!
MODULE mo_jsb_lctlib_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lctlib, lctlib

  TYPE t_lctlib ! Tropical evergreen trees
  !----------------------------------------------------------------------------------------
  ! LANDCOVER LIBRARY
  !
  ! $HeadURL: https://svn.zmaw.de/svn/cosmos/branches/cosmos-landveg/util/running/adjunct_files/jsbach/lctlib_nlct21.def $
  ! $Revision: 7624 $
  !----------------------------------------------------------------------------------------
  !
  ! Contains additional specifications for each landcover type
  !
  ! The data refer to the following landcover classification:
  !
  !  1: Glacier
  !  2: Tropical evergreen trees
  !  3: Tropical deciduous trees
  !  4: Extra-tropical evergreen trees
  !  5: Extra-tropical deciduous trees
  !  6: Temperate broadleaf evergreen trees
  !  7: Temperate broadleaf deciduous trees
  !  8: Coniferous evergreen trees
  !  9: Coniferous deciduous trees
  ! 10: Raingreen shrubs
  ! 11: Deciduous shrubs
  ! 12: C3 grass
  ! 13: C4 grass
  ! 14: Pasture
  ! 15: C3 pasture
  ! 16: C4 pasture
  ! 17: Tundra
  ! 18: Swamp
  ! 19: Crops
  ! 20: C3 crop
  ! 21: C4 crop


    ! R: statt POINTER lieber PARAMETER damit ich call add_jsb_var nicht zur Initialisierung benötige
    ! R: Namen verändert damit Benennung innerhalb von JSBACH4 konsistent ist: MaxLAI zu LaiMax
    REAL(wp) ::                                                   &
           ! LctNumber: Landcover type index numbers (not used in the model)
                                   LctNumber = 2,                 &
           ! LandcoverClass (bare soil: 0; glacier: 1; lake: 2; natural forest: 3; natural grassland: 4; other natural vegetation: 5; crops: 6; pastures: 7)
                                   LandcoverClass = 3,            &
           ! PhenologyType: Phenology type (none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grass: 4; crop: 5)
                                   !PhenologyType = 0,             &
                                   Phenotype = 3,                 &
           ! NitrogenScalingFlag: Is nitrogen scaling needed? (no: 0; yes: 1)
                                   NitrogenScalingFlag = 1,       &
           ! C4flag: Photosynthetic pathway (C3: 0; C4: 1)
                                   C4flag = 0,                    &
!     LOGICAL, POINTER, DIMENSION(:)    ::                       &
!                                         C4flag                   ! Photosynthetic pathway, habe ich auf logical umgestellt.
!                                                                  ! R: könnte auch Real (C3: 0; C4: 1) sein damit der Photosyn. Path
!                                                                  !    erweiterbar ist. Z.B. wenn man Tiles einführt die nicht eindeutig
!                                                                  !    zuzurodnen sind und man daher eine Mischung von C3 und C4 einbauen möchte.
!
           ! CarboxRate: Maximum carboxylation rate at 25 Celsius [1.E-6 * Mol(CO2)/m^2/s] (s. Kattge et al. 2009)
                                   CarboxRate = 39.0,             & ! R: sollte in carbox_rate_max_leaf umbenannt werden
           ! ETransport: Maximum electron transport rate at 25 Celsius [1.E-6 * Mol/m^2/s] (Jmax=1.9*Vmax for C3 plants)
                                   ETransport = 74.1,             &
           ! VegHeight: Vegetation height [m]; for non-vegetation value is arbitrary
                                   VegHeight = 30.0,              &
           ! VegRoughness: Vegetation roughness length [m]; for non-vegetation value is arbitrary
           ! but should be larger than zero to avoid division by zero
                                   VegRoughness = 2.0,            &
           ! MinVegRoughness: Vegetation roughness length at LAI=1 [m]; for non-vegetation value is arbitrary
           ! but should be larger than zero to avoid division by zero
                                   MinVegRoughness = 2.0,         &
           ! MaxVegRoughness: Vegetation roughness length at LAI=LAI_max [m]; for non-vegetation value is arbitrary
           ! but should be larger than zero to avoid division by zero
                                   MaxVegRoughness = 5.0,         &
           ! Frac_npp_2_woodPool: Maximum fraction of NPP used for buildup of woody parts
                                   Frac_npp_2_woodPool = 0.3,     &
           ! Frac_npp_2_reservePool: Maximum fraction of NPP put into the storage pool
                                   Frac_npp_2_reservePool = 0.05, &
           ! Frac_npp_2_exudates: Maximum fraction of NPP put into the root exudates
                                   Frac_npp_2_exudates = 0.05,    &
           ! Frac_green_2_herbivory: Maximum fraction of green cpool put into grazing
                                   Frac_green_2_herbivory = 0.000822,  &
           ! Frac_C_litter_green2atmos: Fraction of Carbon from heterotrophic respiration that is emitted to the atmosphere.
           ! (The rest enters slow pool.)
                                   Frac_C_litter_green2atmos = 0.85, &
           ! Tau_Cpool_litter_leaf: Life time of leaf litter [days]
                                   Tau_Cpool_litter_leaf = 870.0,    &
           ! Tau_Cpool_litter_wood: Life time of woody litter [days]
                                   Tau_Cpool_litter_wood = 10950.0,  &
           ! Tau_Cpool_woods: Time scale of the Cpool_woods and vegetation dynamics [years]
                                   Tau_Cpool_woods = 30.0,           &
           ! LAI_shed_constant: Time in which leaves are constantly shedded [days-1]
                                   LAI_shed_constant = 0.000342,     &
           ! Max_C_content_woods: Maximum carbon content in woody parts [mol(C)/m^2]
           ! (values from IPCC-TAR = , S. 192)
                                   Max_C_content_woods = 2997.25,    &
           ! reserveC2leafC: Factor by which the maximum carbon content of the reserve pool is larger than MAXLAI /SLA
                                   reserveC2leafC = 2.0,             &
           ! frac_wood_2_onsite: Fraction of wood that goes into the anthropogenically controlled onSite woody pool
           ! on converting natural PFTs to agricultural PFTs
                                   frac_wood_2_onsite = 0.6,         &
           ! frac_wood_2_paper: Fraction of wood that goes into the paper (intermediately longlived) pool
           ! on converting natural PFTs to agricultural PFTs
                                   frac_wood_2_paper = 0.4,          &
           ! frac_wood_2_construction: Fraction of wood that goes into the construction
           ! (e.g. furniture = , houses = longlived) pool on converting natural PFTs to agricultural PFTs
                                   frac_wood_2_construction = 0.0,   &
           ! LaiMax: Maximum LAI used in the LogoP scheme [m2/m2]
                                   LaiMax = 7.0,                     &
           ! StemArea: Area of stem and branches [m2/m2] (this is an optical parameter for the snow masking by forests)
                                   StemArea = 1.5,                   &
           ! specificLeafArea_C: Carbon content per leaf area in [m2(leaf)/mol(Carbon)]
                                   specificLeafArea_C = 0.264,       &
           ! ClumpinessFactor: Factor to calculate vegetation clumpiness:
           ! veg_ratio=veg_ratio_max*(1-exp(-LAI_max/ClumpinessFactor)); arbitrary for non-vegetation but needs to be >0 to avoid division by zero
                                   ClumpinessFactor = 2.0,           &
           ! AlbedoCanopyVIS: Albedo of the canopy in the visible range
                                   AlbedoCanopyVIS = 0.03,           &
           ! AlbedoCanopyNIR: Albedo of the canopy in the near infrared range
                                   AlbedoCanopyNIR = 0.22,           &
           ! AlbedoLitterVIS: Albedo of litter in the visible range
                                   AlbedoLitterVIS = 0.09,           &
           ! AlbedoLitterNIR: Albedo of litter in the near infrared range
                                   AlbedolitterNIR = 0.16,           &
           ! AlbedoSnowVisMin: Minimum snow albedo in the visible range
                                   AlbedoSnowVisMin = 0.52,       &
           ! AlbedoSnowVisMax: Maximum snow albedo in the visible range
                                   AlbedoSnowVisMax = 0.9,        &
           ! AlbedoSnowNirMin: Minimum snow albedo in the NIR range
                                   AlbedoSnowNirMin = 0.3,        &
           ! AlbedoSnowNirMax: Maximum snow albedo in the NIR range
                                   AlbedoSnowNirMax = 0.65,       &
           ! Albedo of Snow: only used with the old echam albedo scheme (USE_ALBEDO=.FALSE.)
                                   AlbedoSnowMin = 0.4,          &
                                   AlbedoSnowMax = 0.8,          &
           ! DYNAMIC_PFT: Flag to indicate those PFTs that are subject to dynamics
                                   DYNAMIC_PFT = 1,               &
           ! WOODY_PFT: Flag to indicate woody type PFTs (in contrast to grasses)
                                   WOODY_PFT = 1,                 &
           ! PASTURE_PFT: Flag to indicate pasture (in contrast to crops)
                                   PASTURE_PFT = 0,               &
           ! BCLIMIT_MIN_COLD_mmTEMP: PFT-specific minimum coldest monthly mean temperature
                                   BCLIMIT_MIN_COLD_mmTEMP = 15.5, &
           ! BCLIMIT_MAX_COLD_mmTEMP: PFT-specific maximum coldest monthly mean temperature
                                   BCLIMIT_MAX_COLD_mmTEMP = 1000.0,  &
           !  BCLIMIT_MAX_WARM_mmTEMP: PFT-specific maximum warmest monthly mean temperature
                                   BCLIMIT_MAX_WARM_mmTEMP = 1000.0,  &
           !  BCLIMIT_MIN_TEMPRANGE: PFT-specific 20-year average min warmest - coldest month temperature range
                                   BCLIMIT_MIN_TEMPRANGE = -1000.0,   &
           ! BCLIMIT_MIN_GDD: PFT-specific minimum growing degree days (at or above 5 deg C)
                                   BCLIMIT_MIN_GDD = 0.0,         &
           ! GDD_BASE: PFT-specific GDD base
                                   GDD_BASE = 5.0,                &
           ! UPPER_TLIM: PFT-specific upper limit of warmest-month temperature (used to calculate gdd_upper_tlim)
                                   UPPER_TLIM = 0.0,              &
           !
           !  Parameters for Knorr Phenology
           !
           ! KNORR_TAU_W: Time before leaf shedding [days]
                                   KNORR_TAU_W = 300.0,             &
           ! KNORR_T_PHI: Temperature trigger for leaf growth [deg C]
                                   KNORR_T_PHI = 0.00,            &
           ! KNORR_T_R: Spread (sigma) of T_phi [deg C]
                                   KNORR_T_R = 0.00,              &
           ! KNORR_DAY_C: Day-length at leaf shedding [hours]
                                   KNORR_DAY_C = 0.00,            &
           ! KNORR_DAY_R: Spread (sigma) of Day_c [hours]
                                   KNORR_DAY_R = 0.00,            &
           ! KNORR_K_L: Inverse of leaf longevity [days-1]
                                   KNORR_K_L = 0.07,              &
           ! KNORR_LEAF_GROWTH_RATE: Initial leaf growth rate [days-1]
                                   KNORR_LEAF_GROWTH_RATE = 0.37, &
           ! KNORR_MAX_LAI: Maximum LAI
                                   KNORR_MAX_LAI = 4.2,           &
           !
           !  Parameters for Yasso soil carbon
           !
           ! LitVeg_coef: coefficient to distribute litter into 5 classes of chemical composition
                                   LitVeg_coef = 0.6500,          &
           ! Sum (has to be 1) 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
           ! LeafLit_coef: coefficient to distribute leaf litter into 5 classes of chemical composition
                                   LeafLit_coef = 0.4651,         &
           ! Sum (has to be 1) 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
           ! WoodLit_coef: coefficient to distribute woody litter into 5 classes of chemical composition
                                   WoodLit_coef = 0.6500,         &
           ! Sum (has to be 1) 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
           ! WoodLitterSize:
                                   WoodLitterSize = 4,           &
           !
           !  Parameters for Spitfire fire algorithm (Thonicke et al. 2010)
           !
           ! moist_extinction: moisture of extinction used to determine rate of spread in spitfire units: [-]
                                   moist_extinction = 0.2,         &
           ! fuel_dens: bulk fuel density [kg/m3]
                                   fuel_dens = 25,                &
           ! f parameter to compute the flame length
                                   flame_length_f = 0.1487,      &
           !  crown length parameters
                                   crown_length = 0.33,          &
           !  bark parameter 1
                                   bark_par1 = 0.0301,           &
           !  bark parameter 2
                                   bark_par2 = 0.0281,           &
           !  RCK
                                   rck = 1.00,                   &
           !  probability for mortality
                                   mort_prob = 3.00

  END TYPE t_lctlib ! Tropical evergreen trees


!   TYPE t_lctlib ! Glacier!
!   !----------------------------------------------------------------------------------------
!   ! LANDCOVER LIBRARY
!   !
!   ! $HeadURL: https://svn.zmaw.de/svn/cosmos/branches/cosmos-landveg/util/running/adjunct_files/jsbach/lctlib_nlct21.def $
!   ! $Revision: 7624 $
!   !----------------------------------------------------------------------------------------
!   !
!   ! Contains additional specifications for each landcover type
!   !
!   ! The data refer to the following landcover classification:
!   !
!   !  1: Glacier
!   !  2: Tropical evergreen trees
!   !  3: Tropical deciduous trees
!   !  4: Extra-tropical evergreen trees
!   !  5: Extra-tropical deciduous trees
!   !  6: Temperate broadleaf evergreen trees
!   !  7: Temperate broadleaf deciduous trees
!   !  8: Coniferous evergreen trees
!   !  9: Coniferous deciduous trees
!   ! 10: Raingreen shrubs
!   ! 11: Deciduous shrubs
!   ! 12: C3 grass
!   ! 13: C4 grass
!   ! 14: Pasture
!   ! 15: C3 pasture
!   ! 16: C4 pasture
!   ! 17: Tundra
!   ! 18: Swamp
!   ! 19: Crops
!   ! 20: C3 crop
!   ! 21: C4 crop
!
!
!
!
!
!     ! R: statt POINTER lieber PARAMETER damit ich call add_jsb_var nicht zur Initialisierung benötige
!     REAL, PARAMETER, DIMENSION(:) ::                              &
!            ! LctNumber: Landcover type index numbers (not used in the model)
!                                    LctNumber = 1,                 &
!            ! LandcoverClass (bare soil: 0; glacier: 1; lake: 2; natural forest: 3; natural grassland: 4; other natural vegetation: 5; crops: 6; pastures: 7)
!                                    LandcoverClass = 1,            &
!            ! PhenologyType: Phenology type (none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grass: 4; crop: 5)
!                                    !PhenologyType = 0,             &
!                                    Phenotype = 0,                 &
!            ! NitrogenScalingFlag: Is nitrogen scaling needed? (no: 0; yes: 1)
!                                    NitrogenScalingFlag = 0,       &
!            ! C4flag: Photosynthetic pathway (C3: 0; C4: 1)
!                                    C4flag = 0,                    &
! !     LOGICAL, POINTER, DIMENSION(:)    ::                       &
! !                                         C4flag                   ! Photosynthetic pathway, habe ich auf logical umgestellt.
! !                                                                  ! R: könnte auch Real (C3: 0; C4: 1) sein damit der Photosyn. Path
! !                                                                  !    erweiterbar ist. Z.B. wenn man Tiles einführt die nicht eindeutig
! !                                                                  !    zuzurodnen sind und man daher eine Mischung von C3 und C4 einbauen möchte.
! !
!            ! CarboxRate: Maximum carboxylation rate at 25 Celsius [1.E-6 * Mol(CO2)/m^2/s] (s. Kattge et al. 2009)
!                                    CarboxRate = 0.0,              &
!            ! ETransport: Maximum electron transport rate at 25 Celsius [1.E-6 * Mol/m^2/s] (Jmax=1.9*Vmax for C3 plants)
!                                    ETransport = 0.0,              &
!            ! VegHeight: Vegetation height [m]; for non-vegetation value is arbitrary
!                                    VegHeight = 0.0,               &
!            ! VegRoughness: Vegetation roughness length [m]; for non-vegetation value is arbitrary
!            ! but should be larger than zero to avoid division by zero
!                                    VegRoughness = 0.005,          &
!            ! MinVegRoughness: Vegetation roughness length at LAI=1 [m]; for non-vegetation value is arbitrary
!            ! but should be larger than zero to avoid division by zero
!                                    MinVegRoughness = 0.005,       &
!            ! MaxVegRoughness: Vegetation roughness length at LAI=LAI_max [m]; for non-vegetation value is arbitrary
!            ! but should be larger than zero to avoid division by zero
!                                    MaxVegRoughness = 0.010,       &
!            ! Frac_npp_2_woodPool: Maximum fraction of NPP used for buildup of woody parts
!                                    Frac_npp_2_woodPool = 0.0,     &
!            ! Frac_npp_2_reservePool: Maximum fraction of NPP put into the storage pool
!                                    Frac_npp_2_reservePool = 0.0,  &
!            ! Frac_npp_2_exudates: Maximum fraction of NPP put into the root exudates
!                                    Frac_npp_2_exudates = 0.0,     &
!            ! Frac_green_2_herbivory: Maximum fraction of green cpool put into grazing
!                                    Frac_green_2_herbivory = 0.0,  &
!            ! Frac_C_litter_green2atmos: Fraction of Carbon from heterotrophic respiration that is emitted to the atmosphere.
!            ! (The rest enters slow pool.)
!                                    Frac_C_litter_green2atmos = 1.0, &
!            ! Tau_Cpool_litter_leaf: Life time of leaf litter [days]
!                                    Tau_Cpool_litter_leaf = 1.0,   &
!            ! Tau_Cpool_litter_wood: Life time of woody litter [days]
!                                    Tau_Cpool_litter_wood = 10950.0, &
!            ! Tau_Cpool_woods: Time scale of the Cpool_woods and vegetation dynamics [years]
!                                    Tau_Cpool_woods = 1.0,         &
!            ! LAI_shed_constant: Time in which leaves are constantly shedded [days-1]
!                                    LAI_shed_constant = 0.0,       &
!            ! Max_C_content_woods: Maximum carbon content in woody parts [mol(C)/m^2]
!            ! (values from IPCC-TAR = , S. 192)
!                                    Max_C_content_woods = 0.0,     &
!            ! reserveC2leafC: Factor by which the maximum carbon content of the reserve pool is larger than MAXLAI /SLA
!                                    reserveC2leafC = 2.0,          &
!            ! frac_wood_2_onsite: Fraction of wood that goes into the anthropogenically controlled onSite woody pool
!            ! on converting natural PFTs to agricultural PFTs
!                                    frac_wood_2_onsite = 1.0,      &
!            ! frac_wood_2_paper: Fraction of wood that goes into the paper (intermediately longlived) pool
!            ! on converting natural PFTs to agricultural PFTs
!                                    frac_wood_2_paper = 0.0,       &
!            ! frac_wood_2_construction: Fraction of wood that goes into the construction
!            ! (e.g. furniture = , houses = longlived) pool on converting natural PFTs to agricultural PFTs
!                                    frac_wood_2_construction = 0.0, &
!            ! LaiMax: Maximum LAI used in the LogoP scheme [m2/m2]
!                                    LaiMax = 0.0,                  &
!            ! StemArea: Area of stem and branches [m2/m2] (this is an optical parameter for the snow masking by forests)
!                                    StemArea = 0.0,                &
!            ! specificLeafArea_C: Carbon content per leaf area in [m2(leaf)/mol(Carbon)]
!                                    specificLeafArea_C = 0.0,      &
!            ! ClumpinessFactor: Factor to calculate vegetation clumpiness:
!            ! veg_ratio=veg_ratio_max*(1-exp(-LAI_max/ClumpinessFactor)); arbitrary for non-vegetation but needs to be >0 to avoid division by zero
!                                    ClumpinessFactor = 1.0,        &
!            ! AlbedoCanopyVIS: Albedo of the canopy in the visible range
!                                    AlbedoCanopyVIS = 0.0,         &
!            ! AlbedoCanopyNIR: Albedo of the canopy in the near infrared range
!                                    AlbedoCanopyNIR = 0.0,         &
!            ! AlbedoLitterVIS: Albedo of litter in the visible range
!                                    AlbedoLitterVIS = 0.0,         &
!            ! AlbedoLitterNIR: Albedo of litter in the near infrared range
!                                    AlbedolitterNIR = 0.0,         &
!            ! AlbedoSnowVisMin: Minimum snow albedo in the visible range
!                                    AlbedoSnowVisMin = 0.78,       &
!            ! AlbedoSnowVisMax: Maximum snow albedo in the visible range
!                                    AlbedoSnowVisMax = 0.9,        &
!            ! AlbedoSnowNirMin: Minimum snow albedo in the NIR range
!                                    AlbedoSnowNirMin = 0.4,        &
!            ! AlbedoSnowNirMax: Maximum snow albedo in the NIR range
!                                    AlbedoSnowNirMax = 0.75,       &
!            ! Albedo of Snow: only used with the old echam albedo scheme (USE_ALBEDO=.FALSE.)
!                                    AlbedoSnowMin = 0.75,          &
!                                    AlbedoSnowMax = 0.85,          &
!            ! DYNAMIC_PFT: Flag to indicate those PFTs that are subject to dynamics
!                                    DYNAMIC_PFT = 0,               &
!            ! WOODY_PFT: Flag to indicate woody type PFTs (in contrast to grasses)
!                                    WOODY_PFT = 0,                 &
!            ! PASTURE_PFT: Flag to indicate pasture (in contrast to crops)
!                                    PASTURE_PFT = 0,               &
!            ! BCLIMIT_MIN_COLD_mmTEMP: PFT-specific minimum coldest monthly mean temperature
!                                    BCLIMIT_MIN_COLD_mmTEMP = -1000.0, &
!            ! BCLIMIT_MAX_COLD_mmTEMP: PFT-specific maximum coldest monthly mean temperature
!                                    BCLIMIT_MAX_COLD_mmTEMP = 1000.0,  &
!            !  BCLIMIT_MAX_WARM_mmTEMP: PFT-specific maximum warmest monthly mean temperature
!                                    BCLIMIT_MAX_WARM_mmTEMP = 1000.0,  &
!            !  BCLIMIT_MIN_TEMPRANGE: PFT-specific 20-year average min warmest - coldest month temperature range
!                                    BCLIMIT_MIN_TEMPRANGE = -1000.0,   &
!            ! BCLIMIT_MIN_GDD: PFT-specific minimum growing degree days (at or above 5 deg C)
!                                    BCLIMIT_MIN_GDD = 0.0,         &
!            ! GDD_BASE: PFT-specific GDD base
!                                    GDD_BASE = 5.0,                &
!            ! UPPER_TLIM: PFT-specific upper limit of warmest-month temperature (used to calculate gdd_upper_tlim)
!                                    UPPER_TLIM = 0.0,              &
!            !
!            !  Parameters for Knorr Phenology
!            !
!            ! KNORR_TAU_W: Time before leaf shedding [days]
!                                    KNORR_TAU_W = 0.0,             &
!            ! KNORR_T_PHI: Temperature trigger for leaf growth [deg C]
!                                    KNORR_T_PHI = 0.00,            &
!            ! KNORR_T_R: Spread (sigma) of T_phi [deg C]
!                                    KNORR_T_R = 0.00,              &
!            ! KNORR_DAY_C: Day-length at leaf shedding [hours]
!                                    KNORR_DAY_C = 0.00,            &
!            ! KNORR_DAY_R: Spread (sigma) of Day_c [hours]
!                                    KNORR_DAY_R = 0.00,            &
!            ! KNORR_K_L: Inverse of leaf longevity [days-1]
!                                    KNORR_K_L = 0.00,              &
!            ! KNORR_LEAF_GROWTH_RATE: Initial leaf growth rate [days-1]
!                                    KNORR_LEAF_GROWTH_RATE = 0.37, &
!            ! KNORR_MAX_LAI: Maximum LAI
!                                    KNORR_MAX_LAI = 4.2,           &
!            !
!            !  Parameters for Yasso soil carbon
!            !
!            ! LitVeg_coef: coefficient to distribute litter into 5 classes of chemical composition
!                                    LitVeg_coef = 0.6500,          &
!            ! Sum (has to be 1) 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
!            ! LeafLit_coef: coefficient to distribute leaf litter into 5 classes of chemical composition
!                                    LeafLit_coef = 1.0000,         &
!            ! Sum (has to be 1) 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
!            ! WoodLit_coef: coefficient to distribute woody litter into 5 classes of chemical composition
!                                    WoodLit_coef = 1.0000,         &
!            ! Sum (has to be 1) 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000
!            ! WoodLitterSize:
!                                    WoodLitterSize = 4,           &
!            !
!            !  Parameters for Spitfire fire algorithm (Thonicke et al. 2010)
!            !
!            ! moist_extinction: moisture of extinction used to determine rate of spread in spitfire units: [-]
!                                    moist_extinction = 0,         &
!            ! fuel_dens: bulk fuel density [kg/m3]
!                                    fuel_dens = 0,                &
!            ! f parameter to compute the flame length
!                                    flame_length_f = 0.0000,      &
!            !  crown length parameters
!                                    crown_length = 0.00,          &
!            !  bark parameter 1
!                                    bark_par1 = 0.0000,           &
!            !  bark parameter 2
!                                    bark_par2 = 0.0000,           &
!            !  RCK
!                                    rck = 0.00,                   &
!            !  probability for mortality
!                                    mort_prob = 0.00
!
!   END TYPE t_lctlib  ! Glacier


  TYPE(t_lctlib), SAVE                    :: lctlib

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lct_memory_class'


#endif
END MODULE mo_jsb_lctlib_class
