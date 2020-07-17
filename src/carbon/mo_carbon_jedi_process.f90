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
MODULE mo_carbon_jedi_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_lctlib_class, ONLY: t_lctlib_element

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_Ctradeoff, yasso, get_per_tile,jedi_water_stress_factor,jedi_soil_depth_to_layers_1d!,drop_expensive_tissue
  !PUBLIC :: , relocate_carbon_fire, relocate_carbon_damage
  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_jedi_process'

CONTAINS

  FUNCTION jedi_soil_depth_to_layers_1d(depth, dz) RESULT(depth_l)

    REAL(wp), INTENT(in) :: depth(:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: depth_l(SIZE(depth,1), SIZE(dz))

    REAL(wp) :: bottom, bottom_m1
    INTEGER  :: nsoil, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':jedi_soil_depth_to_layer_1d'

    nsoil = SIZE(dz)

    depth_l(:,:) = 0._wp
    bottom = 0._wp
    bottom_m1 = 0._wp
    DO i=1,nsoil
      bottom = bottom + dz(i)
      WHERE (depth(:) >= bottom)
        depth_l(:,i) = dz(i)
      ELSEWHERE (depth(:) > bottom_m1)
        depth_l(:,i) = depth(:) - bottom_m1
      END WHERE
      bottom_m1 = bottom_m1 + dz(i)
    END DO

  END FUNCTION jedi_soil_depth_to_layers_1d



 ELEMENTAL FUNCTION jedi_water_stress_factor ( &
    & w_soil, w_soil_max, w_soil_wilt_fract) RESULT(water_stress_factor)

    ! TODO
    ! TBD: Maybe change later so that it uses geographically dependend critical values (like wilting point)

    REAL(wp), INTENT(in) :: w_soil            !< Soil water content
    REAL(wp), INTENT(in) :: w_soil_max        !< Soil water content at field capacity
    
    REAL(wp), INTENT(in) :: w_soil_wilt_fract !< Fraction of max. soil water content at wilting point
    REAL(wp)             :: water_stress_factor

    REAL(wp) :: w_crit, w_wilt                       !< Soil water content at critical/wilting point

    
    w_wilt = w_soil_max * w_soil_wilt_fract

    IF (w_soil_max - w_wilt > 0._wp) THEN
      water_stress_factor = MAX (0._wp, &
        &                        MIN (1._wp, (w_soil - w_wilt) / (w_soil_max - w_wilt) ))
    ELSE
      water_stress_factor = 0._wp
    END IF

 END FUNCTION jedi_water_stress_factor



  ! Calculations to allocate NPP to the plant C pools.
  ELEMENTAL SUBROUTINE calc_Ctradeoff( &
    & GPP_yDayMean    ,           & ! in
    ! JeDi plant trait parameters <<<<<<<<<<
    & GrowRspdMoist,              & ! in
    & GrowRspdTemp ,              & ! in
    & CritiGrowTemp,              & ! in
    & EndGrowTemp,              & !in
    & GermFrac     ,              & ! in
    & Allo2Reproduc,              & ! in
    & Allo2Abground,              & ! in
    & Allo2Blground,              & ! in
    & Allo2Storage ,              & ! in
    & RelAllo2Abgrd,              & ! in
    & RelAllo2Blgrd,              & ! in
    & TurnoverStruc,              & ! in
    & TurnoverLeaf ,              & ! in
    & SenesRspdNPP ,              & ! in
    & NitroStatus  ,              & ! in
    & JeDiSLA,                     &  ! in
    & LandcoverClass,          &  !in
    ! JeDi plant trait parameters >>>>>>>>>>
    & iniyear,                    &     ! in
    & currentyear,                & !in
    & time_steps_per_day,         & ! in
    & faPAR_prev_noon,              & !in
    & one_counter,                      & !inout
    & c_storage,                  & ! inout
    & c_leaves,                   & ! inout
    & c_fnroot,                   & ! inout
    & c_csroot,                   & ! inout
    & c_seed,                     & ! inout
    & c_stem,                     & ! inout
    & ff_grow_W,                  & ! in
    & ff_grow_T,                  & ! in
    & ff_grow_G,                  & ! in
    & ff_dTdt,                       & !in
    & Flg_grow,                    & ! out
    & Flg_grow_G ,              & ! out
    & ff_npp,                     & ! inout
    & Flg_die,                    & ! inout
    & Flg_senesc,                 & ! out
    & flx_germ,                   & ! out
    & adj_fraction,                &  !out
    & jedi_light_constraint,         &  !out
    & MaintRes_yDayMean,          & ! inout
    & GrowthResp,                 & ! out
    & AutoRes,                    & ! out
    & CUE,                          & !out carbon use efficiency
    & root_depth,                 & ! out root_depth computed from coarse root carbon
    & C_balan,                    & ! out C_balan = GPP - Mainres - Growthres - litter
    & C_delta,                    & ! out C_delta = New total carbon - old total carbon
    & c_total,                    & ! out
    &  c_Res_N,                    & ! in
    &  c_CN_woody,               & ! in
    &  c_sapwood,                   & ! in
    &  temp_yDayMean,        & ! in
    & c_litter_green_ag,          & ! inout
    & c_litter_green_bg,          & ! inout
    & c_litter_wood_ag,           & ! inout
    & c_litter_wood_bg,           & ! inout
 !
    & soilResp_rate,              & ! out
 !
    & cflux_c_greenwood_2_litter, & ! out
 !
    & NPP_act_yDayMean,           & ! out
 !
    ! Yasso variables
    & temp2_30d,                  & ! in, optional, necessary with yasso
    & precip_30d,                 & ! in, optional, necessary with yasso
    !
    & c_acid_ag1,                 & ! inout, optional, necessary with yasso
    & c_water_ag1,                & ! inout, optional, necessary with yasso
    & c_ethanol_ag1,              & ! inout, optional, necessary with yasso
    & c_nonsoluble_ag1,           & ! inout, optional, necessary with yasso
    & c_acid_bg1,                 & ! inout, optional, necessary with yasso
    & c_water_bg1,                & ! inout, optional, necessary with yasso
    & c_ethanol_bg1,              & ! inout, optional, necessary with yasso
    & c_nonsoluble_bg1,           & ! inout, optional, necessary with yasso
    & c_humus_1,                  & ! inout, optional, necessary with yasso
    & c_acid_ag2,                 & ! inout, optional, necessary with yasso
    & c_water_ag2,                & ! inout, optional, necessary with yasso
    & c_ethanol_ag2,              & ! inout, optional, necessary with yasso
    & c_nonsoluble_ag2,           & ! inout, optional, necessary with yasso
    & c_acid_bg2,                 & ! inout, optional, necessary with yasso
    & c_water_bg2,                & ! inout, optional, necessary with yasso
    & c_ethanol_bg2,              & ! inout, optional, necessary with yasso
    & c_nonsoluble_bg2,           & ! inout, optional, necessary with yasso
    & c_humus_2,                  & ! inout, optional, necessary with yasso
    !
    & LeafLit_coef1,              & ! in, optional, necessary with yasso
    & LeafLit_coef2,              & ! in, optional, necessary with yasso
    & LeafLit_coef3,              & ! in, optional, necessary with yasso
    & LeafLit_coef4,              & ! in, optional, necessary with yasso
    & LeafLit_coef5,              & ! in, optional, necessary with yasso
    & WoodLit_coef1,              & ! in, optional, necessary with yasso
    & WoodLit_coef2,              & ! in, optional, necessary with yasso
    & WoodLit_coef3,              & ! in, optional, necessary with yasso
    & WoodLit_coef4,              & ! in, optional, necessary with yasso
    & WoodLit_coef5,              & ! in, optional, necessary with yasso
    & WoodLitterSize              & ! in, optional, necessary with yasso
    & )

    USE mo_carbon_constants,       ONLY : fract_green_aboveGround, fract_wood_aboveGround
    USE mo_jsb_physical_constants, ONLY : tmelt     ! tmelt =273.15_wp
    USE mo_carbon_constants,       ONLY : sec_per_day, days_per_year, sec_per_year
    USE mo_math_constants,         ONLY: pi
    ! Argument
    
    REAL(wp), INTENT(in)    :: GPP_yDayMean           ! Yesterdays mean GPP rate [mol(C)/m^2 s]
    REAL(wp), INTENT(in)     :: GrowRspdMoist         ! Jedi parameter -- P01 Growth response time to moisture conditions 
    REAL(wp), INTENT(in)     :: GrowRspdTemp          ! Jedi parameter -- P02 Growth response time to temperature conditions 
    REAL(wp), INTENT(in)     :: CritiGrowTemp         ! Jedi parameter -- P03 Critical temperature for growth        
    REAL(wp), INTENT(in)     :: EndGrowTemp        
    REAL(wp), INTENT(in)     :: GermFrac              ! Jedi parameter -- P04 Germination fraction                      
    REAL(wp), INTENT(in)     :: Allo2Reproduc         ! Jedi parameter -- P05 Allocation to reproduction                
    REAL(wp), INTENT(in)     :: Allo2Abground         ! Jedi parameter -- P06 Allocation to above ground growth           
    REAL(wp), INTENT(in)     :: Allo2Blground         ! Jedi parameter -- P07 Allocation to below ground growth           
    REAL(wp), INTENT(in)     :: Allo2Storage          ! Jedi parameter -- P08 Allocation to storage                       
    REAL(wp), INTENT(in)     :: RelAllo2Abgrd         ! Jedi parameter -- P09 Relative allocation to aboveground growth   
    REAL(wp), INTENT(in)     :: RelAllo2Blgrd         ! Jedi parameter -- P10 Relative alloction to belowground growth   
    REAL(wp), INTENT(in)     :: TurnoverStruc         ! Jedi parameter -- P11 Turnover time of structural pools           
    REAL(wp), INTENT(in)     :: TurnoverLeaf          ! Jedi parameter -- P12 Turnover time of leaf and fine root pools   
    REAL(wp), INTENT(in)     :: SenesRspdNPP          ! Jedi parameter -- P13 Senescence response time to productivity conditions 
    REAL(wp), INTENT(in)     :: NitroStatus           ! Jedi parameter -- P15 Plant nitrogen status    
    REAL(wp), INTENT(in)     :: JeDiSLA      
    INTEGER, INTENT(in)     :: LandcoverClass             
    INTEGER,  INTENT(in)     :: iniyear      
    INTEGER,  INTENT(in)     :: currentyear  
    INTEGER,  INTENT(in)     :: time_steps_per_day
    REAL(wp), INTENT(inout) :: one_counter
    !
    ! JeDi carbon allocation 
    REAL(wp), INTENT(inout) :: c_storage            ! Stored carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]                     
    REAL(wp), INTENT(inout) :: c_leaves           ! Leaves carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_seed             ! Seed carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]                                                                                                                                                   
    REAL(wp), INTENT(inout) :: c_fnroot           ! Fine root carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_csroot           ! Coarse root carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]                                                     
    REAL(wp), INTENT(inout) :: c_stem             ! Stem carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]        
    REAL(wp), INTENT(in) :: ff_grow_T,ff_dTdt
    REAL(wp), INTENT(in) :: ff_grow_G
    REAL(wp), INTENT(in) :: ff_grow_W 
    REAL(wp), INTENT(out)   ::  Flg_grow
    REAL(wp), INTENT(out)   ::  Flg_grow_G
    REAL(wp), INTENT(inout) :: ff_npp
    REAL(wp), INTENT(out)   :: flx_germ  
    REAL(wp), INTENT(out)   :: Flg_senesc
    REAL(wp), INTENT(inout) :: Flg_die    
    REAL(wp), INTENT(out)   :: adj_fraction  
    REAL(wp), INTENT(in)    :: faPAR_prev_noon    
    REAL(wp), INTENT(out)   :: jedi_light_constraint                                                  
    !                                                     
    REAL(wp), INTENT(inout) :: c_litter_green_ag  ! Above ground green litter C-pool: updated on output [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_litter_green_bg  ! Below ground green litter C-pool: updated on output [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_litter_wood_ag   ! Above ground woody litter C-pool: updated on output [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_litter_wood_bg   ! Below ground woody litter C-pool: updated on output [mol(C)/m^2(canopy)]
 
    !
    REAL(wp), INTENT(out)   ::  soilResp_rate          ! Soil (=heterotrophic) respiration rate  [mol(C)/m^2 s] Note: this is a loss,
                                                     !    therefore negative!
    REAL(wp), INTENT(out)   ::  cflux_c_greenwood_2_litter
    REAL(wp), INTENT(out)   ::  NPP_act_yDayMean       ! Actual NPP after N-limitation and excess carbon drop [mol(C)/m^2 s]
    REAL(wp), INTENT(inout) ::  MaintRes_yDayMean
    REAL(wp), INTENT(out)   ::  GrowthResp
    REAL(wp),INTENT(IN)      ::   c_Res_N     
    REAL(wp),INTENT(IN)      ::  c_CN_woody
    REAL(wp),INTENT(IN)      ::  c_sapwood  
    REAL(wp),INTENT(IN)      ::   temp_yDayMean     
    REAL(wp), INTENT(out)   ::  AutoRes
    REAL(wp), INTENT(out)   ::  C_balan
    REAL(wp), INTENT(out)   ::  C_delta
    REAL(wp), INTENT(out)   ::  c_total
    REAL(wp), INTENT(out)   ::  CUE
    REAL(wp), INTENT(out)   ::  root_depth
    !
    ! added for thonicke fire algorithm (spitfire) needed to correctly split the carbon pools into the different fuel classes
    !real(wp),optional,intent(inout) :: fract_litter_wood_new   ! new fraction in above ground wood litter pool

    ! Meteorology for Yasso
    REAL(wp), OPTIONAL, INTENT(in)    :: temp2_30d              ! 30 day mean temperature
    REAL(wp), OPTIONAL, INTENT(in)    :: precip_30d             ! 30 day mean precipitation

    ! Yasso pools
    ! Size class 1; green litter
    !   Aboveground
    REAL(wp), OPTIONAL, INTENT(inout) :: c_acid_ag1        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_water_ag1       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_ethanol_ag1     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_nonsoluble_ag1  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground
    REAL(wp), OPTIONAL, INTENT(inout) :: c_acid_bg1        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_water_bg1       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_ethanol_bg1     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_nonsoluble_bg1  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_humus_1         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2; woody litter
    !   Aboveground
    REAL(wp), OPTIONAL, INTENT(inout) :: c_acid_ag2        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_water_ag2       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_ethanol_ag2     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_nonsoluble_ag2  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground
    REAL(wp), OPTIONAL, INTENT(inout) :: c_acid_bg2        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_water_bg2       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_ethanol_bg2     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_nonsoluble_bg2  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), OPTIONAL, INTENT(inout) :: c_humus_2         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Parameters for Yasso
       ! Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(wp), OPTIONAL, INTENT(in)    :: LeafLit_coef1          ! fraction going to the acid soluble pools
    REAL(wp), OPTIONAL, INTENT(in)    :: LeafLit_coef2          ! fraction going to the water soluble pools
    REAL(wp), OPTIONAL, INTENT(in)    :: LeafLit_coef3          ! fraction going to the ethanol soluble pools
    REAL(wp), OPTIONAL, INTENT(in)    :: LeafLit_coef4          ! fraction going to the non soluble pools
    REAL(wp), OPTIONAL, INTENT(in)    :: LeafLit_coef5          ! fraction going to the humus pool
       ! Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(wp), OPTIONAL, INTENT(in)    :: WoodLit_coef1          ! fraction going to the acid soluble pools
    REAL(wp), OPTIONAL, INTENT(in)    :: WoodLit_coef2          ! fraction going to the water soluble pools
    REAL(wp), OPTIONAL, INTENT(in)    :: WoodLit_coef3          ! fraction going to the ethanol soluble pools
    REAL(wp), OPTIONAL, INTENT(in)    :: WoodLit_coef4          ! fraction going to the non soluble pools
    REAL(wp), OPTIONAL, INTENT(in)    :: WoodLit_coef5          ! fraction going to the humus pool
    REAL(wp), OPTIONAL, INTENT(in)    :: WoodLitterSize         ! size of coarse debris



    ! locals
    
    ! JeDi carbon allocation 

    
    !real(wp) :: r8_tiny = 0.0123456
    REAL(wp), PARAMETER  :: PA0_seed = 1.0_wp / (365._wp * 12.01_wp) ! Virtual carbon for seed pool
    
    REAL(wp) :: A_seed , A_storage
    REAL(wp) :: AG_struct, AG_growth
    REAL(wp) :: BG_struct, BG_growth
    REAL(wp) :: Flg_seed
    REAL(wp) :: Flg_grow_T, Flg_grow_W
    REAL(wp) :: tau_stem, tau_leaves, tau_fnroot, tau_csroot, tau_seed
    REAL(wp) :: ff_germ_Bsoil
    REAL(wp) :: T_degC, fract_veg
    
    REAL(wp), PARAMETER  :: c_GRes_leaves    = 1.0_wp / (1.0_wp + 0.33_wp)
    REAL(wp), PARAMETER  :: c_GRes_fnroot    = 1.0_wp / (1.0_wp + 0.33_wp)
    REAL(wp), PARAMETER  :: c_GRes_wood      = 1.0_wp / (1.0_wp + 0.25_wp)
    REAL(wp), PARAMETER  :: c_GRes_csroot    = 1.0_wp / (1.0_wp + 0.49_wp)
    REAL(wp), PARAMETER  :: c_GRes_seed      = 1.0_wp / (1.0_wp + 0.66_wp)    
    REAL(wp), PARAMETER  :: c_tau_seed       = 3._wp * 365._wp !  seed decomposition time scale (years), value is taken from JeDi 
    REAL(wp), PARAMETER  :: c_tau_sen        = 10._wp 
    REAL(wp), PARAMETER  :: c_tau_storage     = 1._wp * 365._wp ! decay of the storage pool (1-year), PHS decided arbitrarily
    REAL(wp), parameter  :: cTiny            = 1.E-8_wp        
    
    REAL(wp), PARAMETER  :: leaf_growth_fract = 0.03_wp   
    REAL(wp), PARAMETER  :: fapar_max= 0.96_wp
    REAL(wp) :: allo_max_leaves,zAlloc_leaves_tmp
    
    REAL(wp) :: MaintRes
    REAL(wp) :: NPP
    REAL(WP) :: f_leafroot,f_finetissues,f_stemcsroot                              
    REAL(WP) :: fract_green_AG_jedi, fract_wood_AG_jedi ! USED FOR YASSOR 
    REAL(WP) :: zAlloc, zAlloc_leaves, zAlloc_stem, zAlloc_seed, zAlloc_fnroot, zAlloc_csroot
    REAL(WP) :: zLoss_leaves, zLoss_stem, zLoss_seed, zLoss_fnroot, zLoss_csroot, zLoss_storage 
    REAL(WP) :: c_total_temp
    REAL(wp) :: C_2_litter_woodPool
    !CBAL REAL(wp) :: d_c_slow
    REAL(wp) :: leaf_shedding 
    REAL(wp) :: alpha_mod
    REAL(wp) :: NPP_2_rootExudates ! PHS : can be delete later
    
     !CBAL REAL(wp) :: Cflx_slow_2_atmos
    REAL(wp) :: Cflx_faeces_2_LG, Cflx_faeces_2_atm
    ! R: only for spitfire
    !REAL(wp) :: litter_wood_new, litter_wood_total

    ! local variables for Yasso
    REAL(wp), DIMENSION(18) :: Yasso_io_pools
    REAL(wp), DIMENSION(2)  :: Weather
    REAL(wp), DIMENSION(16) :: Yasso_out
    REAL(wp), DIMENSION(5)  :: LeafLit_coefV
    REAL(wp), DIMENSION(5)  :: WoodLit_coefV
    
    REAL(wp)                :: leafLitter      ! amount of leafLitter
    REAL(wp)                :: sizeZero
    
    
    REAL(wp) :: decompRate_litter_green

    REAL(wp) :: redFact_Nlimitation
     ! R: Normally parameters should be in mo_carbon_constants, but the following parameters are used only in this subroutine!
    !    Therefore I put them here. Is there any not structural advantage to put them into mo_carbon_constants?
    ! Parameter for Subroutine calc_Cpools
    REAL(wp), PARAMETER :: Q10                    = 1.8_wp          ! Empirical parameter for temperature dependence of
                                                                    !   heterotrophic respiration rate
    REAL(wp), PARAMETER :: kappa                  = 1.0_wp          ! Empirical parameter for soil moisture dependence of
                                                                    !   heterotrophic respiration rate
    REAL(wp), PARAMETER :: referenceTemp_Q10      = Tmelt           ! Reference temperature in the Q10 formula [Kelvin]
    REAL(wp), PARAMETER :: tau_c_reserve          = 1.0_wp *365._wp ! time constant by which reserve pool is depreciated [days]

    REAL(wp), PARAMETER :: tau_c_slow             = 100._wp*365._wp ! time constant by which slow pool is depreciated  [days]

   ! REAL(wp), PARAMETER :: greenC2leafC           = 4.0_wp          ! ratio of carbon green pool (leaves, fine roots, starches,
                                                                    !    sugars) to leaf carbon
    REAL(wp), PARAMETER :: fract_C_litter_wood2atmos = 0.2_wp       ! fraction of heterotrophic loss of woody litter emitted to the
    REAL(wp), PARAMETER :: alpha_critical         = 0.35_wp         ! critical value of relative soil moisture below which
                                                                    !    heterotrophic respiration almost stops
    REAL(wp), PARAMETER :: alpha_min              = 0.1_wp          ! effective lowest value for soil moisture entering the
  
    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    soilResp_rate            = 0.0_wp
    !soilResp_rate_pot        = 0.0_wp
    NPP_act_yDayMean         = 0.0_wp
      
    

    ! Initializations


    cflux_c_greenwood_2_litter = 0.0_wp



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
    leafLitter      = 0.0_wp           ! amount of leafLitter
    sizeZero        = 0.0_wp
    
    flx_germ   = 0._wp
    AutoRes    = 0._wp  
    
    GrowthResp = 0._wp
    zAlloc     = 0._wp
    zLoss_seed = 0._wp
    !r8_tiny = tiny(r8_tiny)  ! Maybe it is different in our station. This is the smallest positive number in the data type of r8_tiny ?/????
    
     
    
    AutoRes = MaintRes_yDayMean * sec_per_day
    NPP     = GPP_yDayMean * sec_per_day  - AutoRes
    
    c_total_temp = c_storage + c_seed + c_leaves + c_stem + c_csroot + c_fnroot 
    
   
   
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Growing condition !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
    
    
    
    Flg_grow_T = 0._wp
    Flg_grow_W = 0._wp
    Flg_grow_G = 0._wp
    
  !  IF ( ff_grow_T     .ge.  CritiGrowTemp   ) Flg_grow_T = 1.0_wp
    
    IF ( ff_grow_T     .ge.  EndGrowTemp   ) THEN
            Flg_grow_T = 1.0_wp
    ELSE IF ( ff_grow_T     .gt.  CritiGrowTemp  .and.  ff_grow_T     .le.  EndGrowTemp  ) THEN
        IF ( ff_dTdt .ge. 0._wp  ) THEN
            Flg_grow_T = 1.0_wp
        ELSE
            Flg_grow_T = 0.0_wp
        END IF
    ELSE 
            Flg_grow_T = 0.0_wp
    END IF
    
    IF ( ff_grow_W     .ge. 0.5_wp           ) Flg_grow_W = 1.0_wp
    IF ( ff_grow_G     .ge. 0.5_wp           ) Flg_grow_G = 1.0_wp  
    
    
    Flg_grow   = Flg_grow_T * Flg_grow_W 
    Flg_grow_G = Flg_grow_T * Flg_grow_G
    
    IF( Flg_grow_G .GT. 0._wp .AND.  one_counter .lt. 1._wp  )  THEN
           
      c_storage    = 100
      one_counter =one_counter +1
    
    ENDIF
   
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Germination condition !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
   
    
    
    IF ( Flg_die .lt. 1._wp )  THEN  ! If the species is alive
    
  
    flx_germ   = Flg_grow_G * GermFrac * c_seed ! if germination is favorable then carbon flows from seed to storage pool 
    c_storage  = c_storage + flx_germ
    c_seed     = c_seed    - flx_germ  
    
    
    zLoss_seed =  c_seed   / c_tau_seed  ! Why the decay of seed is only computed outside the loop?
    c_seed     =  c_seed   - zLoss_seed  !  need to think about when the seed will decay maybe only when its after 150 year of simulation or when its not spinup simulation??
    
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Carbon allocation factors !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
    
    fract_veg = 1.0_wp  ! Dont know what this control, need to check later  ????????
     
  
    
    AG_growth  = Flg_grow * Allo2Abground
    BG_growth  = Flg_grow * Allo2Blground
    
    AG_struct  = Flg_grow * RelAllo2Abgrd * fract_veg
    BG_struct  = Flg_grow * RelAllo2Blgrd * fract_veg
    
 
 
    
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Carbon allocation !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
    
    ! Compute JeDi carbon allocation of different C(tissue) pools
    
    
      IF (   c_storage .gt. 0._wp ) THEN  
	    Flg_seed = 0._wp
    
        IF ( NPP  .gt. 0._wp  ) Flg_seed = 1.0_wp ! OR NPP_yDayMean ???
        
        A_seed    = Flg_grow * Flg_seed * Allo2Reproduc
    
        IF( c_storage .gt. 0._wp )  zLoss_storage = c_storage / c_tau_storage
      
        c_storage = c_storage + NPP - zLoss_storage
        
     
        
        
        
        zAlloc        = MAX(0._wp,c_storage)  !  To make sure that no negative c_storage distributed to the pools
        
        !PHS: 10.10.2019  I introduce a break to stop rapid growth of leaf
        ! I define the initial amount of carbon that is added to the leaf pool to be 4g C.
        ! Afterwards, the maximal allocation to c_leaves can only be 3% c_leaves of the previous day. (leaf_growth_fract)
        ! Likewise, the other tissues' allocation are constrained by the proportional of the c_leaves allocation.  (adj_fraction)
        IF ( c_leaves .EQ. 0._wp )  THEN
           allo_max_leaves = 4._wp / 12.01_wp
        ELSE
           allo_max_leaves = c_leaves * leaf_growth_fract
        END IF
        !PHS: 9.6.2020 FAPAR limitation on leaf growth
           
         IF ( faPAR_prev_noon .LT. fapar_max ) THEN
         !PHS: 15.7.2020 try to  use different light constraint for grass 
          IF         (LandcoverClass .EQ. 3) THEN                    
           jedi_light_constraint = COS(  0.5_wp *  pi   * (   (faPAR_prev_noon / fapar_max  )**30._wp) )
           ELSE IF (LandcoverClass .EQ. 4) THEN
             jedi_light_constraint = COS(  0.5_wp *  pi   * (   (faPAR_prev_noon / 0.83_wp  )**30._wp) )
          END IF   
        ELSE
            jedi_light_constraint = 0._wp
        END IF
        
        zAlloc_leaves_tmp = zAlloc * AG_growth * ( 1._wp - AG_struct ) * c_GRes_leaves  
        zAlloc_leaves = MIN( zAlloc_leaves_tmp, allo_max_leaves )  
        
        IF (zAlloc_leaves .NE. 0._wp  ) THEN
            adj_fraction  =  zAlloc_leaves / zAlloc_leaves_tmp
        ELSE
            adj_fraction  = 1._wp
        END IF
        
        zAlloc_leaves =  zAlloc_leaves  * jedi_light_constraint
        
        
        zAlloc_seed   = zAlloc * A_seed                                          * c_GRes_seed     *   adj_fraction   * jedi_light_constraint
        !zAlloc_leaves = zAlloc * AG_growth * ( 1._wp - AG_struct ) * c_GRes_leaves
        zAlloc_stem   = zAlloc * AG_growth * AG_struct                  * c_GRes_wood     *   adj_fraction   * jedi_light_constraint
        zAlloc_csroot = zAlloc * BG_growth * BG_struct                  * c_GRes_csroot    *   adj_fraction    * jedi_light_constraint
        zAlloc_fnroot = zAlloc * BG_growth * ( 1._wp - BG_struct ) * c_GRes_fnroot     *   adj_fraction    * jedi_light_constraint
        
        c_storage = c_storage - ( A_seed + AG_growth + BG_growth ) * zAlloc   *   adj_fraction   * jedi_light_constraint
        
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Growth respiration !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
        
        
        GrowthResp =  ( zAlloc_leaves  / c_GRes_leaves ) * ( 1.0_wp - c_GRes_leaves ) &
        &                  + ( zAlloc_fnroot  / c_GRes_fnroot ) * ( 1.0_wp - c_GRes_fnroot ) &
        &                  + ( zAlloc_stem   / c_GRes_wood   ) * ( 1.0_wp - c_GRes_wood   ) &
        &                  + ( zAlloc_csroot  / c_GRes_csroot ) * ( 1.0_wp - c_GRes_csroot ) &
        &                  + ( zAlloc_seed    / c_GRes_seed   ) * ( 1.0_wp - c_GRes_seed   ) 
        
        AutoRes = AutoRes + GrowthResp  
        
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Senescence condition !=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
        
        
        
        Flg_senesc = 0._wp
        IF ( (    currentyear .GE. iniyear+20   ) .and. ff_npp .lt. 0._wp .and. NPP .lt. 0._wp) Flg_senesc  = 1.0_wp - Flg_grow   ! See remarks below.
        ! No idea why the flag for senescence is using:  1 - flag_grow  
        ! I took the calculation from JeDi.
        ! I guess senescence happens, when: (1) Non-growing condition (2) NPP<0, which may suggest that plants spend many energy respiring
        !                                   (3) C_storage>0           (4) ff_npp<0 
        
        
        
        tau_stem   = TurnoverStruc
        tau_csroot = TurnoverStruc
        tau_leaves = 1._wp / ( (1._wp / TurnoverLeaf ) + ( 1._wp / TurnoverLeaf ) * Flg_senesc )    ! *             RelSenesAbgd )      !PHS 11.12.2019 play around with the senescence process
        tau_fnroot = 1._wp / ( (1._wp / TurnoverLeaf ) + ( 1._wp / TurnoverLeaf ) * Flg_senesc  )    ! * ( 1._wp -   RelSenesAbgd) )!PHS 11.12.2019 play around with the senescence process
        
        
        
        
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Compute turnover Carbon !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
        
        zLoss_leaves =  c_leaves / tau_leaves
        
        zLoss_stem   =  c_stem   / tau_stem
        zLoss_csroot =  c_csroot  / tau_csroot
        zLoss_fnroot =  c_fnroot   / tau_fnroot
        
               
        
        
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Update Carbon pools !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
        
        
        
        c_seed    = c_seed    + zAlloc_seed        
        c_leaves  = c_leaves  + zAlloc_leaves - zLoss_leaves        
        c_stem    = c_stem   + zAlloc_stem   - zLoss_stem        
        c_csroot  = c_csroot  + zAlloc_csroot - zLoss_csroot       
        c_fnroot  = c_fnroot   + zAlloc_fnroot - zLoss_fnroot
        
        
        
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Death of plants     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
        
       
         
            NPP_act_yDayMean = ( GPP_yDayMean * sec_per_day  - AutoRes ) / (  sec_per_day )
                
            IF         (LandcoverClass .EQ. 3) THEN     
               
               root_depth = sqrt( c_csroot * 2250._wp*12.01_wp ) * 0.001_wp
           
           ELSE IF (LandcoverClass .EQ. 4) THEN
 
                root_depth = sqrt( c_fnroot * 550._wp*12.01_wp ) * 0.001_wp
 
           END IF
      
      
          GrowthResp = GrowthResp /  sec_per_day 
          AutoRes    = AutoRes    /  sec_per_day 
          
    
          IF (  GPP_yDayMean .LT. 1.E-12_wp ) THEN
          
              CUE = 0._wp
              
          ELSE
        
           CUE =NPP_act_yDayMean /SIGN(MAX(ABS(GPP_yDayMean), 1.E-12_wp), GPP_yDayMean)
           CUE = MAX ( 0._wp, CUE )
           
          END IF

        
          IF (  c_storage   .le. cTiny )  then 
              
            zLoss_leaves = zLoss_leaves + c_leaves
            zLoss_stem   = zLoss_stem   + c_stem    + c_storage ! I dont know where the dead storage carbon shall go so I first put it into stem letter pool
            zLoss_csroot = zLoss_csroot + c_csroot  + c_seed ! I dont know where should dead seed carbon go so i put into coarse root's litter pool first
            zLoss_fnroot = zLoss_fnroot + c_fnroot
            
            c_storage = 0._wp   ! Should the storage carbon go to some litter pools?
            c_seed    = 0._wp
            c_leaves  = 0._wp
            c_stem    = 0._wp
            c_csroot  = 0._wp
            c_fnroot  = 0._wp
            Flg_die   = 1._wp
            
            
    
        
        END IF
 
         IF (  Flg_grow  .LT.  1._WP) CUE = 0._wp

        
        
        
          leafLitter          =  zLoss_leaves + zLoss_fnroot + zLoss_storage
          C_2_litter_woodPool =  zLoss_stem   + zLoss_csroot 
          
          
            IF (leafLitter == 0._wp ) THEN 
               fract_green_AG_jedi = 0._wp
            ELSE
               fract_green_AG_jedi =  zLoss_leaves  / leafLitter
            ENDIF
            
            IF (  zLoss_stem   + zLoss_csroot + zLoss_seed == 0._wp ) THEN
               fract_wood_AG_jedi = 0._wp
            ELSE
               fract_wood_AG_jedi  =  zLoss_stem    / ( zLoss_stem   + zLoss_csroot + zLoss_seed)   !need to think about seed litter carbon
            ENDIF
            
          cflux_c_greenwood_2_litter =  leafLitter + zLoss_stem + zLoss_csroot + zLoss_seed 
          
        

        ENDIF
        
       
    
    
    
    
    

    ELSE 
      CUE                          = 0._wp
      GrowthResp                   = 0._wp
      AutoRes                      = 0._wp
      NPP_act_yDayMean             = 0._wp
      leafLitter                   = 0._wp
      C_2_litter_woodPool          = 0._wp
      fract_green_AG_jedi          = 0._wp
      fract_wood_AG_jedi           = 0._wp
      cflux_c_greenwood_2_litter   = 0._wp
      root_depth                   = 0._wp   

    ENDIF

    


   c_total = c_storage + c_leaves + c_stem + c_csroot + c_fnroot 

   C_delta =  ( c_storage + c_seed + c_leaves + c_stem + c_csroot + c_fnroot ) - c_total_temp  ! this is the difference of assimilated carbon added into c-pools in this timestep
   

   C_balan = GPP_yDayMean         *  sec_per_day  &
           & - MaintRes_yDayMean  *  sec_per_day  &
           & - GrowthResp         *  sec_per_day  &
           & - ( zLoss_leaves + zLoss_fnroot + zLoss_stem + zLoss_csroot + zLoss_seed + zLoss_storage) 

   

    !! ============== COMPUTE SOIL RESPIRATION ===================================================================
    ! with yasso:
   ! leafLitter    = cflux_c_greenwood_2_litter - C_2_litter_woodPool + Cflx_faeces_2_LG - Cflx_2_crop_harvest
   
    
    soilResp_rate = 0._wp
    NPP_2_rootExudates = 0._wp
    redFact_Nlimitation = 1._wp
    !! call yasso for leaf litter; this time with nutrient effects
    CALL yasso(Yasso_io_pools(1:9), Weather, leafLitter, LeafLit_coefV,               &
      &        sizeZero, Yasso_out, fract_green_AG_jedi, NPP_2_rootExudates,redFact_Nlimitation)

    c_acid_ag1        = Yasso_out(1)
    c_water_ag1       = Yasso_out(2)
    c_ethanol_ag1     = Yasso_out(3)
    c_nonsoluble_ag1  = Yasso_out(4)
    c_acid_bg1        = Yasso_out(5)
    c_water_bg1       = Yasso_out(6)
    c_ethanol_bg1     = Yasso_out(7)
    c_nonsoluble_bg1  = Yasso_out(8)
    c_humus_1         = Yasso_out(9)

    ! rescale respiration (from leaf litter decomposition) from "per day" to "per second"
    soilResp_rate = Yasso_out(10) / sec_per_day

    !! call yasso for woody litter; this time with nutrient effects
    CALL yasso(Yasso_io_pools(10:18), Weather, C_2_litter_woodPool, WoodLit_coefV,       &
      & WoodLitterSize, Yasso_out, fract_wood_AG_jedi, 0.0_wp, redFact_Nlimitation) ! no exudates

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

   
   
   
   
    !! ============== Litter and herbivory diagnostics (flux from day-1 to s-1)
    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter /sec_per_day
   


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

  END SUBROUTINE calc_Ctradeoff




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
    & cflux_dist_greenreserve_2_soil,    & ! inout, because relocate_carbon_fire is
                                           ! called AFTER relocate_carbon_damage!
    & cflux_dist_woods_2_soil,           & ! inout, because relocate_carbon_fire is
                                           ! called AFTER relocate_carbon_damage!
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
    & cflux_fire_all_2_atm               & ! out
    & )

    USE mo_carbon_constants,  ONLY : fract_green_aboveGround, fract_wood_aboveGround

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
    REAL(wp), INTENT(out)   :: cflux_fire_all_2_atm(:)      ! Carbon immedieatly released by fire per tile [mol(C)m-2(canopy)]
    REAL(wp), INTENT(inout) :: cflux_dist_greenreserve_2_soil(:) ! Amount of carbon relocated by wind and fire damage
                                                                 ! .. to the green litter pools [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   ! Amount of carbon relocated by wind damage and fire
                                                            ! .. to the wood litter pools [mol(C)/m^2(canopy)]
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
    REAL(wp) :: factor(SIZE(burned_fract)) ! Helper

    !! diagnose amount of carbon released from the green and reserve pool to the green litter pools
    !! diagnose amount of carbon released from the wood pool to the woody litter pools
    cflux_dist_greenreserve_2_soil(:) = cflux_dist_greenreserve_2_soil(:) + &
      & (c_green(:) + c_reserve(:)) * (1._wp - fract_green_aboveGround) * burned_fract(:)
    cflux_dist_woods_2_soil(:)  = cflux_dist_woods_2_soil(:) + &
      & c_woods(:) * (1._wp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * burned_fract(:)

    !! determine amount of carbon released from the wood pool, living tissue pools, and litter pools to the atmosphere
    cflux_fire_all_2_atm(:) = ( &
        ! take living parts of plants above ground:
      & (c_green(:) + c_reserve(:)) * fract_green_aboveGround              &
        ! now add the YASSO pools above ground:
      & + c_water_ag1(:) + c_acid_ag1(:) + c_ethanol_ag1(:)                &
      & + c_nonsoluble_ag1(:)                                              &
      & + c_water_ag2(:) + c_acid_ag2(:) + c_ethanol_ag2(:)                &
      & + c_nonsoluble_ag2(:)                                              &
        ! now add the fire affected fraction of wood pools above ground:
      & + c_woods(:) * fire_fract_wood_2_atmos * fract_wood_aboveGround    &
      & ) * burned_fract(:)

    ! Lower down density of aboveground litter pools
    factor(:) = 1._wp - burned_fract(:)
    c_acid_ag1      (:) = factor(:) * c_acid_ag1       (:)
    c_water_ag1     (:) = factor(:) * c_water_ag1      (:)
    c_ethanol_ag1   (:) = factor(:) * c_ethanol_ag1    (:)
    c_nonsoluble_ag1(:) = factor(:) * c_nonsoluble_ag1 (:)

    c_acid_ag2      (:) = factor(:) * c_acid_ag2       (:)
    c_water_ag2     (:) = factor(:) * c_water_ag2      (:)
    c_ethanol_ag2   (:) = factor(:) * c_ethanol_ag2    (:)
    c_nonsoluble_ag2(:) = factor(:) * c_nonsoluble_ag2 (:)

    ! Add green and reserve vegetation carbon to belowground green litter pools
    factor(:) = (1._wp-fract_green_aboveGround) * burned_fract(:) * (c_green(:) + c_reserve(:))
    c_acid_bg1      (:) = c_acid_bg1      (:) + factor(:) * LeafLit_coef(1)
    c_water_bg1     (:) = c_water_bg1     (:) + factor(:) * LeafLit_coef(2)
    c_ethanol_bg1   (:) = c_ethanol_bg1   (:) + factor(:) * LeafLit_coef(3)
    c_nonsoluble_bg1(:) = c_nonsoluble_bg1(:) + factor(:) * LeafLit_coef(4)
    c_humus_1       (:) = c_humus_1       (:) + factor(:) * LeafLit_coef(5)

    ! Add the remaings of burned wood to the aboveground wood pools:
    factor(:) = (1._wp-fire_fract_wood_2_atmos) * fract_wood_aboveGround * burned_fract(:) * c_woods(:)
    c_acid_ag2      (:) = c_acid_ag2      (:) + factor(:) * WoodLit_coef(1)
    c_water_ag2     (:) = c_water_ag2     (:) + factor(:) * WoodLit_coef(2)
    c_ethanol_ag2   (:) = c_ethanol_ag2   (:) + factor(:) * WoodLit_coef(3)
    c_nonsoluble_ag2(:) = c_nonsoluble_ag2(:) + factor(:) * WoodLit_coef(4)
    c_humus_2       (:) = c_humus_2       (:) + factor(:) * WoodLit_coef(5)

    ! Add the remaing of burned wood to the belowground wood pools:
    factor(:) = (1._wp-fract_wood_aboveGround) * burned_fract(:) * c_woods(:)
    c_acid_bg2      (:) = c_acid_bg2      (:) + factor(:) * WoodLit_coef(1)
    c_water_bg2     (:) = c_water_bg2     (:) + factor(:) * WoodLit_coef(2)
    c_ethanol_bg2   (:) = c_ethanol_bg2   (:) + factor(:) * WoodLit_coef(3)
    c_nonsoluble_bg2(:) = c_nonsoluble_bg2(:) + factor(:) * WoodLit_coef(4)
    c_humus_2       (:) = c_humus_2       (:) + factor(:) * WoodLit_coef(5)

    ! Lower down the carbon density of living plant pools
    factor(:) = 1._wp - burned_fract(:)
    c_green(:)   = c_green(:)   * factor(:)
    c_woods(:)   = c_woods(:)   * factor(:)
    c_reserve(:) = c_reserve(:) * factor(:)

  end subroutine relocate_carbon_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (c_green_st, c_reserve_st and
  ! c_woods_st) is partly relocated into the litter pools.
  !
  subroutine relocate_carbon_damage(  &
    & damaged_fract,                  & ! in
    & LeafLit_coef,                   & ! in
    & WoodLit_coef,                   & ! in
    & c_green,                        & ! inout, optional
    & c_woods,                        & ! inout, optional
    & c_reserve,                      & ! inout, optional
    & cflux_dist_greenreserve_2_soil, & ! inout
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

    USE mo_carbon_constants,  ONLY : fract_green_aboveGround, fract_wood_aboveGround

    ! Arguments
    real(wp),intent(in)    :: damaged_fract(:)              ! Fraction of the vegetated area of each tile damaged till the last
                                                                ! call of this routine
    real(wp),intent(in)    :: LeafLit_coef(:)         ! Factor to spread non woody litter into yasso pools [ ]
    real(wp),intent(in)    :: WoodLit_coef(:)         ! Factor to spread woody litterfall into yasso pools [ ]
    real(wp),intent(inout) :: c_green(:)               ! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_woods(:)               ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_reserve(:)             ! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(wp),intent(inout)   :: cflux_dist_greenreserve_2_soil(:) ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
    real(wp),intent(inout)   :: cflux_dist_woods_2_soil(:)        ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
                                                                  ! .. to the green litter pools [mol(C)/m^2(canopy)]

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
    REAL(wp) :: factor(SIZE(damaged_fract))      ! Helper

    ! Initializations

    !! transfer carbon from the wood pool to the wood litter pools nd
    !! transfer carbon from the green and reserve pool to the green litter pool

    cflux_dist_woods_2_soil(:) = cflux_dist_woods_2_soil(:) + c_woods(:) * damaged_fract(:)
    cflux_dist_greenreserve_2_soil(:) = cflux_dist_greenreserve_2_soil(:) + (c_green(:) + c_reserve(:)) * damaged_fract(:)

    ! Add aboveground damaged green and reserve
    factor(:) = (c_green(:) + c_reserve(:)) * fract_green_aboveGround * damaged_fract(:)
    c_acid_ag1      (:) = c_acid_ag1      (:) + factor(:) * LeafLit_coef(1)
    c_water_ag1     (:) = c_water_ag1     (:) + factor(:) * LeafLit_coef(2)
    c_ethanol_ag1   (:) = c_ethanol_ag1   (:) + factor(:) * LeafLit_coef(3)
    c_nonsoluble_ag1(:) = c_nonsoluble_ag1(:) + factor(:) * LeafLit_coef(4)

    ! Add belowground damaged green and reserve
    factor(:) = (c_green(:) + c_reserve(:)) * (1._wp - fract_green_aboveGround) * damaged_fract(:)
    c_acid_bg1      (:) = c_acid_bg1      (:) + factor(:) * LeafLit_coef(1)
    c_water_bg1     (:) = c_water_bg1     (:) + factor(:) * LeafLit_coef(2)
    c_ethanol_bg1   (:) = c_ethanol_bg1   (:) + factor(:) * LeafLit_coef(3)
    c_nonsoluble_bg1(:) = c_nonsoluble_bg1(:) + factor(:) * LeafLit_coef(4)
    factor(:) = (c_green(:) + c_reserve(:)) * damaged_fract(:)
    c_humus_1       (:) = c_humus_1       (:) + factor(:) * LeafLit_coef(5)

    ! Add aboveground damaged wood
    factor(:) = fract_wood_aboveGround  * damaged_fract(:) * c_woods(:)
    c_acid_ag2      (:) = c_acid_ag2      (:) + factor(:) * WoodLit_coef(1)
    c_water_ag2     (:) = c_water_ag2     (:) + factor(:) * WoodLit_coef(2)
    c_ethanol_ag2   (:) = c_ethanol_ag2   (:) + factor(:) * WoodLit_coef(3)
    c_nonsoluble_ag2(:) = c_nonsoluble_ag2(:) + factor(:) * WoodLit_coef(4)

    ! Add belowground damaged wood
    factor(:) = (1._wp - fract_wood_aboveGround)  * damaged_fract(:) * c_woods(:)
    c_acid_bg2      (:) = c_acid_bg2      (:) + factor(:) * WoodLit_coef(1)
    c_water_bg2     (:) = c_water_bg2     (:) + factor(:) * WoodLit_coef(2)
    c_ethanol_bg2   (:) = c_ethanol_bg2   (:) + factor(:) * WoodLit_coef(3)
    c_nonsoluble_bg2(:) = c_nonsoluble_bg2(:) + factor(:) * WoodLit_coef(4)
    factor(:) = c_woods(:) * damaged_fract(:)
    c_humus_2       (:) = c_humus_2       (:) + factor(:) * WoodLit_coef(5)

    ! Lower down the carbon density of living plant pools
    factor(:) = 1._wp - damaged_fract(:)
    c_green(:)   = c_green(:)   * factor(:)
    c_woods(:)   = c_woods(:)   * factor(:)
    c_reserve(:) = c_reserve(:) * factor(:)

  end subroutine relocate_carbon_damage



  PURE SUBROUTINE yasso (Yasso_io_pools, Weather, litter, Lit_coefV,WoodLitterSize, Yasso_out, fract_aboveground, &
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
    REAL(wp),                INTENT(IN)  :: fract_aboveground   ! parameter to seperate litter influx into above- and belowground
                                                               !     part                                                   [ ]
    REAL(wp),                INTENT(IN)  :: NPP_2_rootExudates ! root exudates [daily]                                [mol(C)/m2/day]

    REAL(wp), DIMENSION(16), INTENT(OUT) :: Yasso_out          ! updated Yasso pools(1:9) [mol(c)/m2] & respiration(10) [mol(c)/m2/d]
                                                               ! & Cflx_2_humus (11) [mol(c)/m2/d] ... see update_Cpools for more
                                                               ! information

    REAL(wp), OPTIONAL, INTENT(IN)       ::  redFact_Nlimit

    !local
    ! Yasso pools   [mol(c)/m2]
    REAL(wp)     :: c_acid_ag    ! Aboveground
    REAL(wp)     :: c_water_ag
    REAL(wp)     :: c_ethanol_ag
    REAL(wp)     :: c_nonsoluble_ag
    REAL(wp)     :: c_acid_bg    ! Belowground
    REAL(wp)     :: c_water_bg
    REAL(wp)     :: c_ethanol_bg
    REAL(wp)     :: c_nonsoluble_bg

    REAL(wp)     :: c_humus

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
    REAL(wp)     :: fract_litter_2_acid        ! Fraction of incoming material which enters the acid-soluble fraction of litter
    REAL(wp)     :: fract_litter_2_water       ! Fraction of incoming material which enters the water-soluble fraction of litter
    REAL(wp)     :: fract_litter_2_ethanol     ! Fraction of incoming material which enters the ethanol-soluble fraction of litter
    REAL(wp)     :: fract_litter_2_nonsoluble  ! Fraction of incoming material which enters the nonsoluble fraction of litter
    REAL(wp)     :: fract_litter_2_humus       ! Fraction of incoming material which enters the humus pool

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
    fract_litter_2_acid       = Lit_coefV(1)
    fract_litter_2_water      = Lit_coefV(2)
    fract_litter_2_ethanol    = Lit_coefV(3)
    fract_litter_2_nonsoluble = Lit_coefV(4)
    fract_litter_2_humus      = Lit_coefV(5)

    ! initialize yasso pools (AWEN + humus)
    c_acid_ag        = Yasso_io_pools(1)
    c_water_ag       = Yasso_io_pools(2)
    c_ethanol_ag     = Yasso_io_pools(3)
    c_nonsoluble_ag  = Yasso_io_pools(4)

    c_acid_bg        = Yasso_io_pools(5)
    c_water_bg       = Yasso_io_pools(6)
    c_ethanol_bg     = Yasso_io_pools(7)
    c_nonsoluble_bg  = Yasso_io_pools(8)

    c_humus          = Yasso_io_pools(9)

    ! get the sum of the AWEN pools before decomposition
    Pseudo_litter_green = c_acid_ag          &
                        + c_water_ag         &
                        + c_ethanol_ag       &
                        + c_nonsoluble_ag    &
                        + c_acid_bg          &
                        + c_water_bg         &
                        + c_ethanol_bg       &
                        + c_nonsoluble_bg

    ! Change units of climatic forcing variables
    precip = Weather(2)*sec_per_year/1000._wp  ! mm/s -> m/a
    temp   = Weather(1) - Tmelt                ! K -> C

    ! Calculate annual litter influxes [mol(c)/m2/a]
    Cflx_litter_2_acid       = fract_litter_2_acid        * litter * days_per_year
    Cflx_litter_2_water      = fract_litter_2_water       * litter * days_per_year
    Cflx_litter_2_ethanol    = fract_litter_2_ethanol     * litter * days_per_year
    Cflx_litter_2_nonsoluble = fract_litter_2_nonsoluble  * litter * days_per_year
    Cflx_litter_2_humus      = fract_litter_2_humus       * litter * days_per_year


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
     d_acid       =  redFactor * ref_decomp_rate_acid       * d_temp * d_precip * d_size
     d_water      =  redFactor * ref_decomp_rate_water      * d_temp * d_precip * d_size
     d_ethanol    =  redFactor * ref_decomp_rate_ethanol    * d_temp * d_precip * d_size
     d_nonsoluble =  redFactor * ref_decomp_rate_nonsoluble * d_temp * d_precip * d_size

     ! the decomposition of humus is not limited by nutrients; thus no redFactor
     d_humus      =              ref_decomp_rate_humus      * d_temp * d_precip        ! no size effect on humus

    ! 2. Calculate fluxes

    ! loss fluxes (negative values):
    Cflx_from_acid_ag         = c_acid_ag       * d_acid
    Cflx_from_water_ag        = c_water_ag      * d_water
    Cflx_from_ethanol_ag      = c_ethanol_ag    * d_ethanol
    Cflx_from_nonsoluble_ag   = c_nonsoluble_ag * d_nonsoluble

    Cflx_from_acid_bg         = c_acid_bg       * d_acid
    Cflx_from_water_bg        = c_water_bg      * d_water
    Cflx_from_ethanol_bg      = c_ethanol_bg    * d_ethanol
    Cflx_from_nonsoluble_bg   = c_nonsoluble_bg * d_nonsoluble

    Cflx_from_humus           = c_humus         * d_humus

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
    c_acid_ag           = MAX(0.0_wp,                                         &
                             c_acid_ag                                        & ! old pool
                             + (Cflx_from_acid_ag                             & ! incoming flux from AWEN pools
                               + Cflx_2_acid_ag                               & ! outgoing flux
                               + fract_aboveground * Cflx_litter_2_acid) *dt)   ! incoming flux from litter
    c_water_ag          = MAX(0.0_wp,                                         & ! and so on .....
                             c_water_ag                                       &
                             + (Cflx_from_water_ag                            &
                               + Cflx_2_water_ag                              &
                               + fract_aboveground * Cflx_litter_2_water) *dt)
    c_ethanol_ag        = MAX(0.0_wp,                                         &
                             c_ethanol_ag                                     &
                             + (Cflx_from_ethanol_ag                          &
                               + Cflx_2_ethanol_ag                            &
                               + fract_aboveground * Cflx_litter_2_ethanol) *dt)
    c_nonsoluble_ag     = MAX(0.0_wp,                                         &
                             c_nonsoluble_ag                                  &
                             + (Cflx_from_nonsoluble_ag                       &
                               + Cflx_2_nonsoluble_ag                         &
                               + fract_aboveground * Cflx_litter_2_nonsoluble) *dt)

    c_acid_bg           = MAX(0.0_wp,                                         &
                             c_acid_bg                                        & ! old pool
                             + (Cflx_from_acid_bg                             & ! incoming flux from AWEN pools
                               + Cflx_2_acid_bg                               & ! outgoing flux
                               + (1._wp - fract_aboveground) *                & !
                                 Cflx_litter_2_acid) *dt)                       ! incoming flux from litter
    c_water_bg          = MAX(0.0_wp,                                         & ! ...
                             c_water_bg                                       &
                             + (Cflx_from_water_bg                            &
                               + Cflx_2_water_bg                              &
                               + (1._wp - fract_aboveground) *                &
                                 Cflx_litter_2_water) *dt                     &
                             + NPP_2_rootExudates)                               ! exudates are carbohydrates only and belowground
    c_ethanol_bg        = MAX(0.0_wp,                                         &
                             c_ethanol_bg                                     &
                             + (Cflx_from_ethanol_bg                          &
                               + Cflx_2_ethanol_bg                            &
                               + (1._wp - fract_aboveground) *                &
                                    Cflx_litter_2_ethanol) *dt)
    c_nonsoluble_bg     = MAX(0.0_wp,                                         &
                             c_nonsoluble_bg                                  &
                             + (Cflx_from_nonsoluble_bg                       &
                               + Cflx_2_nonsoluble_bg                         &
                               + (1._wp - fract_aboveground) *                &
                                 Cflx_litter_2_nonsoluble) *dt)

    c_humus          = MAX(0.0_wp,                            &
                             c_humus                          &
                             + (Cflx_from_humus               &
                               + Cflx_2_humus                 &
                               + Cflx_litter_2_humus) *dt)

    ! compute d_litter_green; this is given by the change in AWEN pools
    !          Flx                 Pool(t)                   Pool(t+1)
    !    d_P = ----     &   Flx = ----------   --->   d_P = -----------  - 1
    !          Pool(t)             Pool(t+1)                 Pool(t)

    IF (Pseudo_litter_green .GT. 0.0_wp) THEN
       d_litter_green = (c_acid_ag            &
                         + c_water_ag         &
                         + c_ethanol_ag       &
                         + c_nonsoluble_ag    &
                         + c_acid_bg          &
                         + c_water_bg         &
                         + c_ethanol_bg       &
                         + c_nonsoluble_bg)   &
                         / Pseudo_litter_green - 1._wp
    ELSE
      d_litter_green = 0.0_wp
    ENDIF

    ! 4.1 Write pools into the output array
    Yasso_out(1)           = c_acid_ag
    Yasso_out(2)           = c_water_ag
    Yasso_out(3)           = c_ethanol_ag
    Yasso_out(4)           = c_nonsoluble_ag
    Yasso_out(5)           = c_acid_bg
    Yasso_out(6)           = c_water_bg
    Yasso_out(7)           = c_ethanol_bg
    Yasso_out(8)           = c_nonsoluble_bg
    Yasso_out(9)           = c_humus
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

  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2
  !
  ELEMENTAL SUBROUTINE get_per_tile(  &
    & Cpool_per_tile,                 & ! out
    & Cpool_per_canopy,               & ! in
    & veg_fract_correction,           & ! in
    & fract_fpc_max)                    ! in

    REAL(wp), INTENT(IN) ::  fract_fpc_max, veg_fract_correction
    REAL(wp), INTENT(IN) ::  Cpool_per_canopy
    REAL(wp), INTENT(OUT)::  Cpool_per_tile

    Cpool_per_tile = Cpool_per_canopy * veg_fract_correction * fract_fpc_max

  END SUBROUTINE get_per_tile

#endif
END MODULE mo_carbon_jedi_process
