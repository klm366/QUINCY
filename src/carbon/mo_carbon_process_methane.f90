!>
!! @brief Contains carbon process routines specific for the methane producing soil
!!
!! @author
!!   Stiig Wilkenskjeld based on JSBACH3 code by Thomas Kleinen (2023-03-27)
!!
!! @par Revision History
!!   2023-03-27 : First JSBACH4 version, Stiig Wilkenskjeld, MPI-M

!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
MODULE mo_carbon_process_methane
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_carbon_shared,    ONLY: yasso
  USE mo_carbon_constants, ONLY: molarMassC_kg ! This is a bad place for this constant!
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_methane_soil_c

CONTAINS

  PURE ELEMENTAL SUBROUTINE calc_methane_soil_c( &
    ! General inputs
    &  is_peat, dt,                              &
    ! Input coefficients (fixed)                 
    & LeafLit_coef_acid,                         &
    & LeafLit_coef_water,                        &
    & LeafLit_coef_ethanol,                      &
    & LeafLit_coef_nonsoluble,                   &
    & LeafLit_coef_humus,                        &
    & WoodLit_coef_acid,                         &
    & WoodLit_coef_water,                        &
    & WoodLit_coef_ethanol,                      &
    & WoodLit_coef_nonsoluble,                   &
    & WoodLit_coef_humus,                        &
    & WoodLitterSize,                            &
    & frac_inun_min,                             &
    ! Variable inputs                            
    & temp,                                      &
    & precip,                                    &
    & frac_inun,                                 &
    & wt_peat,                                   &
    & woodLitter,                                &
    & leafLitter,                                &
    & rootExudates,                              &
    ! Yasso pools. Used on input and updated on output
    & c_acid_ag1,                                &
    & c_water_ag1,                               &
    & c_ethanol_ag1,                             &
    & c_nonsoluble_ag1,                          &
    & c_acid_bg1,                                &
    & c_water_bg1,                               &
    & c_ethanol_bg1,                             &
    & c_nonsoluble_bg1,                          &
    & c_humus_1,                                 &
    & c_acid_ag2,                                &
    & c_water_ag2,                               &
    & c_ethanol_ag2,                             &
    & c_nonsoluble_ag2,                          &
    & c_acid_bg2,                                &
    & c_water_bg2,                               &
    & c_ethanol_bg2,                             &
    & c_nonsoluble_bg2,                          &
    & c_humus_2,                                 &
    ! Output fluxes                              
    & c_decomp_humus_1_sum,                      &
    & c_decomp_humus_2_sum,                      &
    & c_into_humus_1_sum,                        &
    & c_into_humus_2_sum,                        &
    & soilResp_rate,                             &
    & soilResp_oxic,                             &
    & soilResp_anoxic,                           &
    & acrotelm_height_ag, acrotelm_height_bg,    &
    & redFact_Nlim                               &
    & )
    USE mo_carbon_constants,  ONLY: fract_green_aboveGround, fract_wood_aboveGround, &
      &                             rho_YCpool_acro, rho_YCpool_cato, Cfrac_YCpool_acro, Cfrac_YCpool_cato, &
      &                             i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble,   &
      &                             i_lctlib_humus
    LOGICAL,            INTENT(in)  :: is_peat                      ! tile is peatland tile
    REAL(wp),           INTENT(in)  :: dt                           ! Carbon time step [s]
    !
    ! Fixed input parameters for Yasso
    !   Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(wp),           INTENT(in)    :: LeafLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(wp),           INTENT(in)    :: LeafLit_coef_water         ! fraction going to the water soluble pools
    REAL(wp),           INTENT(in)    :: LeafLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(wp),           INTENT(in)    :: LeafLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(wp),           INTENT(in)    :: LeafLit_coef_humus         ! fraction going to the humus pool
    !   Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(wp),           INTENT(in)    :: WoodLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(wp),           INTENT(in)    :: WoodLit_coef_water         ! fraction going to the water soluble pools
    REAL(wp),           INTENT(in)    :: WoodLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(wp),           INTENT(in)    :: WoodLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(wp),           INTENT(in)    :: WoodLit_coef_humus         ! fraction going to the humus pool
    REAL(wp),           INTENT(in)    :: WoodLitterSize             ! size of coarse debris
    !   Other constants
    REAL(wp),           INTENT(in)    :: frac_inun_min              ! Minimum inundated fraction to consider wet cell parts

    ! Carbon inputs [mol m-2 s-1]:
    REAL(wp),           INTENT(in)    :: leafLitter
    REAL(wp),           INTENT(in)    :: woodLitter
    REAL(wp),           INTENT(in)    :: rootExudates

    ! Input meteorology for Yasso
    REAL(wp),           INTENT(in)    :: temp        ! multi-day mean temperature [K]
    REAL(wp),           INTENT(in)    :: precip      ! multi-mean precipitation [mm s-1]

    ! Input parameters from topmodel (or other wetland scheme)
    REAL(wp),           INTENT(in)    :: frac_inun   ! Fraction of tile inundated []
    REAL(wp),           INTENT(in)    :: wt_peat     ! Water table on peatlands [m]

    ! Yasso pools [mol(C) m-2(canopy)], input and updated on output
    ! Size class 1; green litter
    !   Aboveground
    REAL(wp),           INTENT(inout) :: c_acid_ag1            ! Acid soluble litter pool
    REAL(wp),           INTENT(inout) :: c_water_ag1           ! Water soluble litter pool
    REAL(wp),           INTENT(inout) :: c_ethanol_ag1         ! Ethanol soluble litter pool
    REAL(wp),           INTENT(inout) :: c_nonsoluble_ag1      ! Non-soluble litter pool
    !   Belowground                                  
    REAL(wp),           INTENT(inout) :: c_acid_bg1            ! Acid soluble litter pool
    REAL(wp),           INTENT(inout) :: c_water_bg1           ! Water soluble litter pool
    REAL(wp),           INTENT(inout) :: c_ethanol_bg1         ! Ethanol soluble litter pool
    REAL(wp),           INTENT(inout) :: c_nonsoluble_bg1      ! Non-soluble litter pool
    REAL(wp),           INTENT(inout) :: c_humus_1             ! Humus litter pool
    ! Size class 2; woody litter                     
    !   Aboveground                                  
    REAL(wp),           INTENT(inout) :: c_acid_ag2            ! Acid soluble litter pool
    REAL(wp),           INTENT(inout) :: c_water_ag2           ! Water soluble litter pool
    REAL(wp),           INTENT(inout) :: c_ethanol_ag2         ! Ethanol soluble litter pool
    REAL(wp),           INTENT(inout) :: c_nonsoluble_ag2      ! Non-soluble litter pool
    !   Belowground                                  
    REAL(wp),           INTENT(inout) :: c_acid_bg2            ! Acid soluble litter pool
    REAL(wp),           INTENT(inout) :: c_water_bg2           ! Water soluble litter pool
    REAL(wp),           INTENT(inout) :: c_ethanol_bg2         ! Ethanol soluble litter pool
    REAL(wp),           INTENT(inout) :: c_nonsoluble_bg2      ! Non-soluble litter pool
    REAL(wp),           INTENT(inout) :: c_humus_2             ! Humus litter pool

    ! Decomposistion fluxes
    REAL(wp),           INTENT(inout) :: c_decomp_humus_1_sum  ! annual cflux sum of humus decomposed (leaf)
    REAL(wp),           INTENT(inout) :: c_decomp_humus_2_sum  ! annual cflux sum of humus decomposed (wood)
    REAL(wp),           INTENT(inout) :: c_into_humus_1_sum    ! annual cflux sum into humus (leaf)
    REAL(wp),           INTENT(inout) :: c_into_humus_2_sum    ! annual cflux sum into humus (wood)

    REAL(wp),           INTENT(out)   :: soilResp_rate         ! soil respiration rate, full
    REAL(wp),           INTENT(out)   :: soilResp_oxic         ! soil respiration rate, oxic part
    REAL(wp),           INTENT(out)   :: soilResp_anoxic       ! soil respiration rate, anoxic part

    ! Diagnostics for methane processes
    REAL(wp),           INTENT(out)   :: acrotelm_height_ag    ! 
    REAL(wp),           INTENT(out)   :: acrotelm_height_bg    ! 
    REAL(wp), OPTIONAL, INTENT(in)    :: redFact_Nlim          ! reduction factor through N limitation

    REAL(wp) :: LeafLit_coefV(5)
    REAL(wp) :: WoodLit_coefV(5)
    REAL(wp) :: Yasso_io_pools(18)
    REAL(wp) :: Yasso_out_oxic(18)
    REAL(wp) :: Yasso_out_anoxic(18)
    REAL(wp) :: weather(2)
    REAL(wp) :: frac_inundated_ag
    REAL(wp) :: frac_inundated_bg
    REAL(wp) :: frac_inundated_YPeat(18)
    REAL(wp) :: acrotelm_C_ag
    REAL(wp) :: acrotelm_C_bg
    REAL(wp) :: catotelm_C
    REAL(wp) :: catotelm_height
    REAL(wp) :: red_Fact_N
    INTEGER  :: i

    IF (PRESENT(redFact_Nlim)) THEN
      red_Fact_N = redFact_Nlim
    ELSE
      red_Fact_N = 1.0_wp
    END IF

    soilResp_rate   = 0._wp
    soilResp_oxic   = 0._wp
    soilResp_anoxic = 0._wp

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

    Weather(1) = temp
    Weather(2) = precip

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

    ! distinguish three cases:
    ! 1) normal: no inundation, no peatland
    ! 2) inundated, but no peatland
    ! 3) peatland
    IF (.NOT. is_peat) THEN ! first do non-peatlands

      ! either without inundation or inundated fraction too small (set up in initialisation part)
      IF (frac_inun <= frac_inun_min) THEN

        ! non-inundated area -- should be the same as in unmodified mo_cbal_cpools

        !! call yasso for leaf litter; this time with nutrient effects
        CALL yasso(dt, Yasso_io_pools(1:9), Weather, leafLitter*dt, LeafLit_coefV,                &
                        0._wp,          Yasso_out_oxic, fract_green_aboveGround, rootExudates*dt, &
                        redFact_Nlimit = red_Fact_N)

        c_acid_ag1           =                        Yasso_out_oxic(1)
        c_water_ag1          =                        Yasso_out_oxic(2)
        c_ethanol_ag1        =                        Yasso_out_oxic(3)
        c_nonsoluble_ag1     =                        Yasso_out_oxic(4)
        c_acid_bg1           =                        Yasso_out_oxic(5)
        c_water_bg1          =                        Yasso_out_oxic(6)
        c_ethanol_bg1        =                        Yasso_out_oxic(7)
        c_nonsoluble_bg1     =                        Yasso_out_oxic(8)
        c_humus_1            =                        Yasso_out_oxic(9)
        soilResp_rate        =                        Yasso_out_oxic(10) / dt
        c_decomp_humus_1_sum = c_decomp_humus_1_sum + Yasso_out_oxic(17)
        c_into_humus_1_sum   = c_into_humus_1_sum   + Yasso_out_oxic(18)

        !! call yasso for woody litter; this time with nutrient effects
        CALL yasso(dt, Yasso_io_pools(10:18), Weather, woodLitter*dt, WoodLit_coefV,           &
                        WoodLitterSize, Yasso_out_oxic, fract_wood_aboveGround, 0.0_wp,        &
                        redFact_Nlimit = red_Fact_N) ! no exudates

        c_acid_ag2           =                        Yasso_out_oxic(1)
        c_water_ag2          =                        Yasso_out_oxic(2)
        c_ethanol_ag2        =                        Yasso_out_oxic(3)
        c_nonsoluble_ag2     =                        Yasso_out_oxic(4)
        c_acid_bg2           =                        Yasso_out_oxic(5)
        c_water_bg2          =                        Yasso_out_oxic(6)
        c_ethanol_bg2        =                        Yasso_out_oxic(7)
        c_nonsoluble_bg2     =                        Yasso_out_oxic(8)
        c_humus_2            =                        Yasso_out_oxic(9)
        soilResp_rate        = soilResp_rate        + Yasso_out_oxic(10) / dt
        c_decomp_humus_2_sum = c_decomp_humus_2_sum + Yasso_out_oxic(17)
        c_into_humus_2_sum   = c_into_humus_2_sum   + Yasso_out_oxic(18)
        
        soilResp_oxic = soilResp_rate
        soilResp_anoxic = 0._wp

      ELSE ! (frac_inun <= frac_inun_min)

        ! YASSO for inundated area
        
        ! in inundated areas, soil C oxic / anoxic distinction is horizontal
        ! therefore use full C pools as input to Yasso, but apply inundated/non-inundated 
        ! fraction to output

        !! call yasso for leaf litter
        !! non-inundated fraction first
        CALL yasso(dt, Yasso_io_pools(1:9), Weather, leafLitter*dt, LeafLit_coefV,              &
                        0._wp,           Yasso_out_oxic, fract_green_aboveGround, rootExudates*dt, &
                        redFact_Nlimit = red_Fact_N, is_anoxic = .FALSE.)
                        
        c_acid_ag1           =                        (1.0_wp-frac_inun) * Yasso_out_oxic(1)
        c_water_ag1          =                        (1.0_wp-frac_inun) * Yasso_out_oxic(2)
        c_ethanol_ag1        =                        (1.0_wp-frac_inun) * Yasso_out_oxic(3)
        c_nonsoluble_ag1     =                        (1.0_wp-frac_inun) * Yasso_out_oxic(4)
        c_acid_bg1           =                        (1.0_wp-frac_inun) * Yasso_out_oxic(5)
        c_water_bg1          =                        (1.0_wp-frac_inun) * Yasso_out_oxic(6)
        c_ethanol_bg1        =                        (1.0_wp-frac_inun) * Yasso_out_oxic(7)
        c_nonsoluble_bg1     =                        (1.0_wp-frac_inun) * Yasso_out_oxic(8)
        c_humus_1            =                        (1.0_wp-frac_inun) * Yasso_out_oxic(9)
        soilResp_oxic        =                        (1.0_wp-frac_inun) * Yasso_out_oxic(10) / dt
        c_decomp_humus_1_sum = c_decomp_humus_1_sum + (1.0_wp-frac_inun) * Yasso_out_oxic(17)
        c_into_humus_1_sum   = c_into_humus_1_sum   + (1.0_wp-frac_inun) * Yasso_out_oxic(18)

        !! inundated fraction second
        CALL yasso(dt, Yasso_io_pools(1:9), Weather, leafLitter*dt, LeafLit_coefV,               &
                        0._wp,          Yasso_out_anoxic, fract_green_aboveGround, rootExudates*dt, &
                        redFact_Nlimit = red_Fact_N, is_anoxic = .TRUE.)
                        
        c_acid_ag1           = c_acid_ag1           +         frac_inun  * Yasso_out_anoxic(1)
        c_water_ag1          = c_water_ag1          +         frac_inun  * Yasso_out_anoxic(2)
        c_ethanol_ag1        = c_ethanol_ag1        +         frac_inun  * Yasso_out_anoxic(3)
        c_nonsoluble_ag1     = c_nonsoluble_ag1     +         frac_inun  * Yasso_out_anoxic(4)
        c_acid_bg1           = c_acid_bg1           +         frac_inun  * Yasso_out_anoxic(5)
        c_water_bg1          = c_water_bg1          +         frac_inun  * Yasso_out_anoxic(6)
        c_ethanol_bg1        = c_ethanol_bg1        +         frac_inun  * Yasso_out_anoxic(7)
        c_nonsoluble_bg1     = c_nonsoluble_bg1     +         frac_inun  * Yasso_out_anoxic(8)
        c_humus_1            = c_humus_1            +         frac_inun  * Yasso_out_anoxic(9)
        soilResp_anoxic      =                                frac_inun  * Yasso_out_anoxic(10) / dt
        c_decomp_humus_1_sum = c_decomp_humus_1_sum +         frac_inun  * Yasso_out_anoxic(17)
        c_into_humus_1_sum   = c_into_humus_1_sum   +         frac_inun  * Yasso_out_anoxic(18)

        !! call yasso for woody litter
        !! non-inundated fraction first
        CALL yasso(dt, Yasso_io_pools(10:18), Weather, woodLitter*dt, WoodLit_coefV,    &
                        WoodLitterSize, Yasso_out_oxic, fract_wood_aboveGround, 0.0_wp, &
                        redFact_Nlimit = red_Fact_N, is_anoxic = .FALSE.) ! no exudates
                        
        c_acid_ag2           =                        (1.0_wp-frac_inun) * Yasso_out_oxic(1)
        c_water_ag2          =                        (1.0_wp-frac_inun) * Yasso_out_oxic(2)
        c_ethanol_ag2        =                        (1.0_wp-frac_inun) * Yasso_out_oxic(3)
        c_nonsoluble_ag2     =                        (1.0_wp-frac_inun) * Yasso_out_oxic(4)
        c_acid_bg2           =                        (1.0_wp-frac_inun) * Yasso_out_oxic(5)
        c_water_bg2          =                        (1.0_wp-frac_inun) * Yasso_out_oxic(6)
        c_ethanol_bg2        =                        (1.0_wp-frac_inun) * Yasso_out_oxic(7)
        c_nonsoluble_bg2     =                        (1.0_wp-frac_inun) * Yasso_out_oxic(8)
        c_humus_2            =                        (1.0_wp-frac_inun) * Yasso_out_oxic(9)
        soilResp_oxic        = soilResp_oxic  +       (1.0_wp-frac_inun) * Yasso_out_oxic(10) / dt
        c_decomp_humus_2_sum = c_decomp_humus_2_sum + (1.0_wp-frac_inun) * Yasso_out_oxic(17)
        c_into_humus_2_sum   = c_into_humus_2_sum   + (1.0_wp-frac_inun) * Yasso_out_oxic(18)

        !! inundated fraction second
        CALL yasso(dt, Yasso_io_pools(10:18), Weather, woodLitter*dt, WoodLit_coefV,      &
                        WoodLitterSize, Yasso_out_anoxic, fract_wood_aboveGround, 0.0_wp, &
                        redFact_Nlimit = red_Fact_N, is_anoxic = .TRUE.) ! no exudates
                        
        c_acid_ag2           = c_acid_ag2           +         frac_inun  * Yasso_out_anoxic(1)
        c_water_ag2          = c_water_ag2          +         frac_inun  * Yasso_out_anoxic(2)
        c_ethanol_ag2        = c_ethanol_ag2        +         frac_inun  * Yasso_out_anoxic(3)
        c_nonsoluble_ag2     = c_nonsoluble_ag2     +         frac_inun  * Yasso_out_anoxic(4)
        c_acid_bg2           = c_acid_bg2           +         frac_inun  * Yasso_out_anoxic(5)
        c_water_bg2          = c_water_bg2          +         frac_inun  * Yasso_out_anoxic(6)
        c_ethanol_bg2        = c_ethanol_bg2        +         frac_inun  * Yasso_out_anoxic(7)
        c_nonsoluble_bg2     = c_nonsoluble_bg2     +         frac_inun  * Yasso_out_anoxic(8)
        c_humus_2            = c_humus_2            +         frac_inun  * Yasso_out_anoxic(9)
        soilResp_anoxic      = soilResp_anoxic      +         frac_inun  * Yasso_out_anoxic(10) / dt
        c_decomp_humus_2_sum = c_decomp_humus_2_sum +         frac_inun  * Yasso_out_anoxic(17)
        c_into_humus_2_sum   = c_into_humus_2_sum   +         frac_inun  * Yasso_out_anoxic(18)

        soilResp_rate = soilResp_oxic + soilResp_anoxic
        
      END IF ! with_inundation

    ELSE ! is_peat

      ! YASSO for peatland case
      acrotelm_C_ag      = 0._wp
      acrotelm_C_bg      = 0._wp
      catotelm_C         = 0._wp
      acrotelm_height_ag = 0._wp
      acrotelm_height_bg = 0._wp
      catotelm_height    = 0._wp

      ! in peatland areas, soil C oxic / anoxic distinction is vertical
      ! therefore use distinguish oxic / anoxic C pools as input to Yasso, do not apply inundated/non-inundated 
      ! fraction to output

      acrotelm_C_ag = SUM(Yasso_io_pools(1:4)) + SUM(Yasso_io_pools(10:13))
      acrotelm_C_bg = SUM(Yasso_io_pools(5:8)) + SUM(Yasso_io_pools(14:17))
      catotelm_C    =     Yasso_io_pools(  9)  +     Yasso_io_pools(   18)

      ! acrotelm surface is defined at the top of the litter layer
      ! acrotelm height is determined from C pool size for ag and bg separately
      ! impose minimum acrotelm height to allow division
      acrotelm_height_ag = MAX(1.e-9_wp,(acrotelm_C_ag * molarMassC_kg * 1000._wp) / (rho_YCpool_acro * Cfrac_YCpool_acro))
      acrotelm_height_bg = MAX(1.e-9_wp,(acrotelm_C_bg * molarMassC_kg * 1000._wp) / (rho_YCpool_acro * Cfrac_YCpool_acro))
      catotelm_height    = MAX(1.e-9_wp,(catotelm_C    * molarMassC_kg * 1000._wp) / (rho_YCpool_cato * Cfrac_YCpool_cato))

      ! anoxic fraction depends on wt position
      ! wt above surface -> everything anoxic
      IF (wt_peat > 0) THEN
        frac_inundated_ag = 1._wp
        frac_inundated_bg = 1._wp
      ELSE
        ! wt below surface, below ag acrotelm
        IF (-wt_peat > acrotelm_height_ag) THEN
          frac_inundated_ag = 0._wp
          frac_inundated_bg = MAX(0._wp,MIN(1._wp,1._wp - (-wt_peat - acrotelm_height_ag)/acrotelm_height_bg))
        ELSE
          ! wt below surface, in ag acrotelm
          frac_inundated_ag = MAX(0._wp,MIN(1._wp,1._wp - (-wt_peat)/acrotelm_height_ag))
          frac_inundated_bg = 1._wp
        END IF
      END IF

      frac_inundated_YPeat( 1: 4) = frac_inundated_ag
      frac_inundated_YPeat( 5: 8) = frac_inundated_bg
      frac_inundated_YPeat(10:13) = frac_inundated_ag
      frac_inundated_YPeat(14:17) = frac_inundated_bg
      frac_inundated_YPeat(    9) = 1._wp
      frac_inundated_YPeat(   18) = 1._wp

      !! call yasso for leaf litter
      !! non-inundated fraction first
      CALL yasso(dt, (1._wp-frac_inundated_YPeat(1:9))*Yasso_io_pools(1:9), Weather, leafLitter*dt, LeafLit_coefV, &
                 0._wp, Yasso_out_oxic, fract_green_aboveGround, rootExudates*dt,                                     &
                 redFact_Nlimit = red_Fact_N, is_anoxic = .FALSE., is_peat = .TRUE.)

      ! There the area modification is done in the yasso call, thus not scaled using inun_frac
      c_acid_ag1           =                        Yasso_out_oxic(1)
      c_water_ag1          =                        Yasso_out_oxic(2)
      c_ethanol_ag1        =                        Yasso_out_oxic(3)
      c_nonsoluble_ag1     =                        Yasso_out_oxic(4)
      c_acid_bg1           =                        Yasso_out_oxic(5)
      c_water_bg1          =                        Yasso_out_oxic(6)
      c_ethanol_bg1        =                        Yasso_out_oxic(7)
      c_nonsoluble_bg1     =                        Yasso_out_oxic(8)
      c_humus_1            =                        Yasso_out_oxic(9)
      soilResp_oxic        =                        Yasso_out_oxic(10) / dt
      c_decomp_humus_1_sum = c_decomp_humus_1_sum + Yasso_out_oxic(17)
      c_into_humus_1_sum   = c_into_humus_1_sum   + Yasso_out_oxic(18)

      !! inundated fraction second
      CALL yasso(dt, frac_inundated_YPeat(1:9)*Yasso_io_pools(1:9), Weather, 0._wp, LeafLit_coefV, &
                      0._wp, Yasso_out_anoxic, fract_green_aboveGround, 0._wp,                     &
                      redFact_Nlimit = red_Fact_N, is_anoxic = .TRUE., is_peat = .TRUE.)
                
      c_acid_ag1           = c_acid_ag1           + Yasso_out_anoxic(1)
      c_water_ag1          = c_water_ag1          + Yasso_out_anoxic(2)
      c_ethanol_ag1        = c_ethanol_ag1        + Yasso_out_anoxic(3)
      c_nonsoluble_ag1     = c_nonsoluble_ag1     + Yasso_out_anoxic(4)
      c_acid_bg1           = c_acid_bg1           + Yasso_out_anoxic(5)
      c_water_bg1          = c_water_bg1          + Yasso_out_anoxic(6)
      c_ethanol_bg1        = c_ethanol_bg1        + Yasso_out_anoxic(7)
      c_nonsoluble_bg1     = c_nonsoluble_bg1     + Yasso_out_anoxic(8)
      c_humus_1            = c_humus_1            + Yasso_out_anoxic(9)
      soilResp_anoxic      =                        Yasso_out_anoxic(10) / dt
      c_decomp_humus_1_sum = c_decomp_humus_1_sum + Yasso_out_anoxic(17)
      c_into_humus_1_sum   = c_into_humus_1_sum   + Yasso_out_anoxic(18)

      !! call yasso for woody litter
      !! non-inundated fraction first
      CALL yasso(dt,(1._wp-frac_inundated_YPeat(10:18))*Yasso_io_pools(10:18), Weather, woodLitter*dt, WoodLit_coefV, &
                 WoodLitterSize, Yasso_out_oxic, fract_wood_aboveGround, 0.0_wp,                                      &
                 redFact_Nlimit = red_Fact_N, is_anoxic = .FALSE., is_peat = .TRUE.) ! no exudates
                
      c_acid_ag2           =                        Yasso_out_oxic(1)
      c_water_ag2          =                        Yasso_out_oxic(2)
      c_ethanol_ag2        =                        Yasso_out_oxic(3)
      c_nonsoluble_ag2     =                        Yasso_out_oxic(4)
      c_acid_bg2           =                        Yasso_out_oxic(5)
      c_water_bg2          =                        Yasso_out_oxic(6)
      c_ethanol_bg2        =                        Yasso_out_oxic(7)
      c_nonsoluble_bg2     =                        Yasso_out_oxic(8)
      c_humus_2            =                        Yasso_out_oxic(9)
      soilResp_oxic        = soilResp_oxic        + Yasso_out_oxic(10) / dt
      c_decomp_humus_2_sum = c_decomp_humus_2_sum + Yasso_out_oxic(17)
      c_into_humus_2_sum   = c_into_humus_2_sum   + Yasso_out_oxic(18)

      !! inundated fraction second
      CALL yasso(dt,frac_inundated_YPeat(10:18)*Yasso_io_pools(10:18), Weather, 0._wp, WoodLit_coefV, &
                 WoodLitterSize, Yasso_out_anoxic, fract_wood_aboveGround, 0.0_wp,                    &
                 redFact_Nlimit = red_Fact_N, is_anoxic = .TRUE., is_peat = .TRUE.) ! no exudates
                
      c_acid_ag2           = c_acid_ag2           + Yasso_out_anoxic(1)
      c_water_ag2          = c_water_ag2          + Yasso_out_anoxic(2)
      c_ethanol_ag2        = c_ethanol_ag2        + Yasso_out_anoxic(3)
      c_nonsoluble_ag2     = c_nonsoluble_ag2     + Yasso_out_anoxic(4)
      c_acid_bg2           = c_acid_bg2           + Yasso_out_anoxic(5)
      c_water_bg2          = c_water_bg2          + Yasso_out_anoxic(6)
      c_ethanol_bg2        = c_ethanol_bg2        + Yasso_out_anoxic(7)
      c_nonsoluble_bg2     = c_nonsoluble_bg2     + Yasso_out_anoxic(8)
      c_humus_2            = c_humus_2            + Yasso_out_anoxic(9)
      soilResp_anoxic      = soilResp_anoxic      + Yasso_out_anoxic(10) / dt
      c_decomp_humus_2_sum = c_decomp_humus_2_sum + Yasso_out_anoxic(17)
      c_into_humus_2_sum   = c_into_humus_2_sum   + Yasso_out_anoxic(18)

      soilResp_rate = soilResp_oxic + soilResp_anoxic

    END IF  ! is_peat

  END SUBROUTINE calc_methane_soil_c
 
#endif
END MODULE mo_carbon_process_methane
