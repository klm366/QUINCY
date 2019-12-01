MODULE mo_carbon_process
 USE mo_kind , ONLY: wp
 PUBLIC :: calc_cpools
 PUBLIC :: relocate_carbon_fire
 PUBLIC :: relocate_carbon_damage
 PUBLIC :: yasso
 PUBLIC :: get_per_tile
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_carbon_process"

CONTAINS
 ELEMENTAL SUBROUTINE calc_cpools ( lai , lai_previous , max_lai ,&
  npp_pot_ydaymean , topsoiltemp , alpha_ydaymean , fract_npp_2_woodpool ,&
  fract_npp_2_reservepool , fract_npp_2_exudates , fract_green_2_herbivory ,&
  tau_c_litter_green , tau_c_litter_wood , tau_c_woods , lai_shed_constant ,&
  fract_c_litter_green2atmos , max_c_content_woods , specific_leaf_area_c ,&
  reservec2leafc , pheno_type , c_green , c_woods , c_reserve , c_crop_harvest&
  , c_litter_green_ag , c_litter_green_bg , c_litter_wood_ag ,&
  c_litter_wood_bg , soilresp_rate , npp_flux_correction , excess_npp ,&
  root_exudates , cflux_c_greenwood_2_litter , cflx_herbivory ,&
  cflux_herb_2_littergreen , cflux_herb_2_atm , npp_act_ydaymean , temp2_30d ,&
  precip_30d , c_acid_ag1 , c_water_ag1 , c_ethanol_ag1 , c_nonsoluble_ag1 ,&
  c_acid_bg1 , c_water_bg1 , c_ethanol_bg1 , c_nonsoluble_bg1 , c_humus_1 ,&
  c_acid_ag2 , c_water_ag2 , c_ethanol_ag2 , c_nonsoluble_ag2 , c_acid_bg2 ,&
  c_water_bg2 , c_ethanol_bg2 , c_nonsoluble_bg2 , c_humus_2 , leaflit_coef1 ,&
  leaflit_coef2 , leaflit_coef3 , leaflit_coef4 , leaflit_coef5 ,&
  woodlit_coef1 , woodlit_coef2 , woodlit_coef3 , woodlit_coef4 ,&
  woodlit_coef5 , woodlittersize )
  USE mo_carbon_constants , ONLY: fract_green_aboveground ,&
   fract_wood_aboveground
  USE mo_jsb_physical_constants , ONLY: tmelt
  USE mo_carbon_constants , ONLY: sec_per_day

  REAL ( KIND= wp ) , INTENT(IN) :: lai
  REAL ( KIND= wp ) , INTENT(IN) :: lai_previous
  REAL ( KIND= wp ) , INTENT(IN) :: max_lai
  REAL ( KIND= wp ) , INTENT(IN) :: npp_pot_ydaymean
  REAL ( KIND= wp ) , INTENT(IN) :: topsoiltemp
  REAL ( KIND= wp ) , INTENT(IN) :: alpha_ydaymean
  REAL ( KIND= wp ) , INTENT(IN) :: fract_npp_2_woodpool
  REAL ( KIND= wp ) , INTENT(IN) :: fract_npp_2_reservepool
  REAL ( KIND= wp ) , INTENT(IN) :: fract_npp_2_exudates
  REAL ( KIND= wp ) , INTENT(IN) :: fract_green_2_herbivory
  REAL ( KIND= wp ) , INTENT(IN) :: tau_c_litter_green
  REAL ( KIND= wp ) , INTENT(IN) :: tau_c_litter_wood
  REAL ( KIND= wp ) , INTENT(IN) :: tau_c_woods
  REAL ( KIND= wp ) , INTENT(IN) :: lai_shed_constant
  REAL ( KIND= wp ) , INTENT(IN) :: fract_c_litter_green2atmos
  REAL ( KIND= wp ) , INTENT(IN) :: max_c_content_woods
  REAL ( KIND= wp ) , INTENT(IN) :: specific_leaf_area_c
  REAL ( KIND= wp ) , INTENT(IN) :: reservec2leafc
  INTEGER , INTENT(IN) :: pheno_type
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_green
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_woods
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_reserve
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_crop_harvest
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_litter_green_ag
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_litter_green_bg
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_litter_wood_ag
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_litter_wood_bg
  REAL ( KIND= wp ) , INTENT(OUT) :: soilresp_rate
  REAL ( KIND= wp ) , INTENT(OUT) :: npp_flux_correction
  REAL ( KIND= wp ) , INTENT(OUT) :: excess_npp
  REAL ( KIND= wp ) , INTENT(OUT) :: root_exudates
  REAL ( KIND= wp ) , INTENT(OUT) :: cflux_c_greenwood_2_litter
  REAL ( KIND= wp ) , INTENT(OUT) :: cflx_herbivory
  REAL ( KIND= wp ) , INTENT(OUT) :: cflux_herb_2_littergreen
  REAL ( KIND= wp ) , INTENT(OUT) :: cflux_herb_2_atm
  REAL ( KIND= wp ) , INTENT(OUT) :: npp_act_ydaymean
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: temp2_30d
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: precip_30d
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_acid_ag1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_water_ag1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_ethanol_ag1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_nonsoluble_ag1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_acid_bg1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_water_bg1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_ethanol_bg1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_nonsoluble_bg1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_humus_1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_acid_ag2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_water_ag2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_ethanol_ag2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_nonsoluble_ag2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_acid_bg2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_water_bg2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_ethanol_bg2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_nonsoluble_bg2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: c_humus_2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: leaflit_coef1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: leaflit_coef2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: leaflit_coef3
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: leaflit_coef4
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: leaflit_coef5
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: woodlit_coef1
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: woodlit_coef2
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: woodlit_coef3
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: woodlit_coef4
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: woodlit_coef5
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: woodlittersize
  LOGICAL :: is_crop
  REAL ( KIND= wp ) :: d_c_litter_green_ag
  REAL ( KIND= wp ) :: d_c_litter_green_bg
  REAL ( KIND= wp ) :: d_c_litter_wood_ag
  REAL ( KIND= wp ) :: d_c_litter_wood_bg
  REAL ( KIND= wp ) :: c_green_max
  REAL ( KIND= wp ) :: c_reserve_optimal
  REAL ( KIND= wp ) :: npp_2_greenpool
  REAL ( KIND= wp ) :: npp_2_woodpool
  REAL ( KIND= wp ) :: npp_2_reservepool
  REAL ( KIND= wp ) :: npp_2_rootexudates
  REAL ( KIND= wp ) :: green_2_herbivory
  REAL ( KIND= wp ) :: c_2_litter_greenpool
  REAL ( KIND= wp ) :: c_2_litter_woodpool
  REAL ( KIND= wp ) :: excess_carbon
  REAL ( KIND= wp ) :: leaf_shedding
  REAL ( KIND= wp ) :: alpha_mod
  REAL ( KIND= wp ) :: c_green_pot
  REAL ( KIND= wp ) :: c_woods_pot
  REAL ( KIND= wp ) :: c_reserve_pot
  REAL ( KIND= wp ) :: c_litter_green_ag_pot
  REAL ( KIND= wp ) :: c_litter_green_bg_pot
  REAL ( KIND= wp ) :: c_litter_wood_ag_pot
  REAL ( KIND= wp ) :: c_litter_wood_bg_pot
  REAL ( KIND= wp ) :: cflx_npp_2_green_pot
  REAL ( KIND= wp ) :: cflx_npp_2_wood_pot
  REAL ( KIND= wp ) :: cflx_npp_2_reserve_pot
  REAL ( KIND= wp ) :: cflx_littergreen_ag_2_atmos_pot
  REAL ( KIND= wp ) :: cflx_littergreen_bg_2_atmos_pot
  REAL ( KIND= wp ) :: cflx_littergreen_ag_2_slow_pot
  REAL ( KIND= wp ) :: cflx_littergreen_bg_2_slow_pot
  REAL ( KIND= wp ) :: cflx_litter_wood_ag_2_atmos_pot
  REAL ( KIND= wp ) :: cflx_litter_wood_bg_2_atmos_pot
  REAL ( KIND= wp ) :: cflx_litter_wood_ag_2_slow_pot
  REAL ( KIND= wp ) :: cflx_litter_wood_bg_2_slow_pot
  REAL ( KIND= wp ) :: cflx_2_crop_harvest
  REAL ( KIND= wp ) :: cflx_crop_harvest_2_atmos
  REAL ( KIND= wp ) :: cflx_faeces_2_lg
  REAL ( KIND= wp ) :: cflx_faeces_2_atm
  REAL ( KIND= wp ) :: yasso_io_pools ( 1 : 18 )
  REAL ( KIND= wp ) :: weather ( 1 : 2 )
  REAL ( KIND= wp ) :: yasso_out ( 1 : 16 )
  REAL ( KIND= wp ) :: leaflit_coefv ( 1 : 5 )
  REAL ( KIND= wp ) :: woodlit_coefv ( 1 : 5 )
  REAL ( KIND= wp ) :: leaflitter
  REAL ( KIND= wp ) :: sizezero
  REAL ( KIND= wp ) :: decomprate_litter_green
  REAL ( KIND= wp ) :: npp_2_green_via_mobilen
  REAL ( KIND= wp ) :: npp_2_wood_via_mobilen
  REAL ( KIND= wp ) :: redfact_nlimitation
  REAL ( KIND= wp ) , PARAMETER :: q10 = 1.8_wp
  REAL ( KIND= wp ) , PARAMETER :: kappa = 1.0_wp
  REAL ( KIND= wp ) , PARAMETER :: referencetemp_q10 = tmelt
  REAL ( KIND= wp ) , PARAMETER :: tau_c_reserve = 1.0_wp * 365._wp
  REAL ( KIND= wp ) , PARAMETER :: tau_c_slow = 100._wp * 365._wp
  REAL ( KIND= wp ) , PARAMETER :: tau_c_crop_harvest = 1._wp * 365._wp
  REAL ( KIND= wp ) , PARAMETER :: greenc2leafc = 4.0_wp
  REAL ( KIND= wp ) , PARAMETER :: fract_c_litter_wood2atmos = 0.2_wp
  REAL ( KIND= wp ) , PARAMETER :: fract_c_faeces2_lg = 0.3_wp
  REAL ( KIND= wp ) , PARAMETER :: fract_c_crop_harvest = 0.5_wp
  REAL ( KIND= wp ) , PARAMETER :: alpha_critical = 0.35_wp
  REAL ( KIND= wp ) , PARAMETER :: alpha_min = 0.1_wp
  REAL ( KIND= wp ) , PARAMETER :: n2o_rate_nitrification = 0.001_wp
  REAL ( KIND= wp ) , PARAMETER :: n2o_rate_denitrification = 0.00125_wp
  REAL ( KIND= wp ) , PARAMETER :: n2o_rate = 0.005_wp
  REAL ( KIND= wp ) , PARAMETER :: sminn_nh4_fraction = 0.4_wp

  soilresp_rate = 0.0_wp
  npp_act_ydaymean = 0.0_wp
  IF ( pheno_type == 5 ) THEN
   is_crop = .TRUE.
  ELSE
   is_crop = .FALSE.
  END IF
  alpha_mod = max ( alpha_min , ( alpha_ydaymean - alpha_critical ) / ( 1._wp&
   - alpha_critical ) )
  npp_flux_correction = 0.0_wp
  c_green_max = 0.0_wp
  c_reserve_optimal = 0.0_wp
  excess_npp = 0.0_wp
  root_exudates = 0.0_wp
  cflux_c_greenwood_2_litter = 0.0_wp
  cflx_herbivory = 0.0_wp
  cflux_herb_2_littergreen = 0.0_wp
  cflux_herb_2_atm = 0.0_wp
  cflx_2_crop_harvest = 0.0_wp
  yasso_io_pools ( 1 ) = c_acid_ag1
  yasso_io_pools ( 2 ) = c_water_ag1
  yasso_io_pools ( 3 ) = c_ethanol_ag1
  yasso_io_pools ( 4 ) = c_nonsoluble_ag1
  yasso_io_pools ( 5 ) = c_acid_bg1
  yasso_io_pools ( 6 ) = c_water_bg1
  yasso_io_pools ( 7 ) = c_ethanol_bg1
  yasso_io_pools ( 8 ) = c_nonsoluble_bg1
  yasso_io_pools ( 9 ) = c_humus_1
  yasso_io_pools ( 10 ) = c_acid_ag2
  yasso_io_pools ( 11 ) = c_water_ag2
  yasso_io_pools ( 12 ) = c_ethanol_ag2
  yasso_io_pools ( 13 ) = c_nonsoluble_ag2
  yasso_io_pools ( 14 ) = c_acid_bg2
  yasso_io_pools ( 15 ) = c_water_bg2
  yasso_io_pools ( 16 ) = c_ethanol_bg2
  yasso_io_pools ( 17 ) = c_nonsoluble_bg2
  yasso_io_pools ( 18 ) = c_humus_2
  weather ( 1 ) = temp2_30d
  weather ( 2 ) = precip_30d
  leaflit_coefv ( 1 ) = leaflit_coef1
  leaflit_coefv ( 2 ) = leaflit_coef2
  leaflit_coefv ( 3 ) = leaflit_coef3
  leaflit_coefv ( 4 ) = leaflit_coef4
  leaflit_coefv ( 5 ) = leaflit_coef5
  woodlit_coefv ( 1 ) = woodlit_coef1
  woodlit_coefv ( 2 ) = woodlit_coef2
  woodlit_coefv ( 3 ) = woodlit_coef3
  woodlit_coefv ( 4 ) = woodlit_coef4
  woodlit_coefv ( 5 ) = woodlit_coef5
  yasso_out ( 1 : 16 ) = 0.0_wp
  leaflitter = 0.0_wp
  sizezero = 0.0_wp
  c_green_max = greenc2leafc * lai / specific_leaf_area_c
  c_reserve_optimal = reservec2leafc * max_lai / specific_leaf_area_c
  leaf_shedding = max ( lai * lai_shed_constant , lai_previous - lai )
  c_2_litter_greenpool = min ( greenc2leafc * leaf_shedding /&
   specific_leaf_area_c , c_green )
  c_green = c_green - c_2_litter_greenpool
  excess_carbon = max ( 0.0_wp , c_green - c_green_max )
  c_2_litter_greenpool = c_2_litter_greenpool + excess_carbon
  c_green = c_green - excess_carbon
  IF ( is_crop ) THEN
   cflx_2_crop_harvest = fract_green_aboveground * fract_c_crop_harvest *&
    c_2_litter_greenpool
   c_litter_green_ag = c_litter_green_ag + fract_green_aboveground * ( 1._wp -&
    fract_c_crop_harvest ) * c_2_litter_greenpool
  ELSE
   c_litter_green_ag = c_litter_green_ag + fract_green_aboveground *&
    c_2_litter_greenpool
  END IF
  c_litter_green_bg = c_litter_green_bg + ( 1._wp - fract_green_aboveground )&
   * c_2_litter_greenpool
  cflux_c_greenwood_2_litter = c_2_litter_greenpool
  green_2_herbivory = c_green * fract_green_2_herbivory
  c_green = c_green - green_2_herbivory
  cflx_faeces_2_lg = green_2_herbivory * fract_c_faeces2_lg
  cflx_faeces_2_atm = green_2_herbivory * ( 1.0_wp - fract_c_faeces2_lg )
  c_litter_green_ag = c_litter_green_ag + cflx_faeces_2_lg
  c_2_litter_woodpool = c_woods / tau_c_woods
  c_litter_wood_ag = c_litter_wood_ag + fract_wood_aboveground *&
   c_2_litter_woodpool
  c_litter_wood_bg = c_litter_wood_bg + ( 1.0_wp - fract_wood_aboveground ) *&
   c_2_litter_woodpool
  c_woods = c_woods - c_2_litter_woodpool
  cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter +&
   c_2_litter_woodpool
  c_2_litter_greenpool = c_reserve / tau_c_reserve
  IF ( is_crop ) THEN
   cflx_2_crop_harvest = cflx_2_crop_harvest + fract_green_aboveground *&
    fract_c_crop_harvest * c_2_litter_greenpool
   c_litter_green_ag = c_litter_green_ag + fract_green_aboveground * ( 1._wp -&
    fract_c_crop_harvest ) * c_2_litter_greenpool
  ELSE
   c_litter_green_ag = c_litter_green_ag + fract_green_aboveground *&
    c_2_litter_greenpool
  END IF
  c_litter_green_bg = c_litter_green_bg + ( 1._wp - fract_green_aboveground )&
   * c_2_litter_greenpool
  c_reserve = c_reserve - c_2_litter_greenpool
  c_crop_harvest = c_crop_harvest + cflx_2_crop_harvest
  cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter +&
   c_2_litter_greenpool
  c_green_pot = c_green
  c_woods_pot = c_woods
  c_reserve_pot = c_reserve
  IF ( npp_pot_ydaymean <= 0.0_wp ) THEN
   npp_2_reservepool = npp_pot_ydaymean * sec_per_day
   c_reserve = c_reserve + npp_2_reservepool
   IF ( c_reserve < 0.0_wp ) THEN
    npp_2_reservepool = npp_2_reservepool - c_reserve
    npp_flux_correction = ( - c_reserve / sec_per_day )
    c_reserve = 0.0_wp
   END IF
   npp_2_reservepool = 0.0_wp
   npp_2_greenpool = 0.0_wp
   npp_2_woodpool = 0.0_wp
   npp_2_rootexudates = 0.0_wp
   c_reserve_pot = c_reserve
  ELSE
   npp_2_woodpool = fract_npp_2_woodpool * npp_pot_ydaymean * sec_per_day
   npp_2_reservepool = fract_npp_2_reservepool * npp_pot_ydaymean *&
    sec_per_day
   npp_2_rootexudates = fract_npp_2_exudates * npp_pot_ydaymean * sec_per_day
   npp_2_greenpool = ( 1.0_wp - fract_npp_2_woodpool - fract_npp_2_reservepool&
    - fract_npp_2_exudates ) * npp_pot_ydaymean * sec_per_day
   c_woods_pot = c_woods_pot + npp_2_woodpool
   excess_carbon = max ( 0.0_wp , c_woods_pot - max_c_content_woods )
   c_woods_pot = c_woods_pot - excess_carbon
   npp_2_greenpool = npp_2_greenpool + excess_carbon
   npp_2_woodpool = npp_2_woodpool - excess_carbon
   IF ( c_reserve_pot < c_reserve_optimal ) THEN
    c_reserve_pot = c_reserve_pot + npp_2_reservepool
    excess_carbon = max ( 0.0_wp , c_reserve_pot - c_reserve_optimal )
    c_reserve_pot = c_reserve_pot - excess_carbon
    npp_2_greenpool = npp_2_greenpool + excess_carbon
    npp_2_reservepool = npp_2_reservepool - excess_carbon
   ELSE
    npp_2_greenpool = npp_2_greenpool + npp_2_reservepool
    npp_2_reservepool = 0.0_wp
   END IF
   c_green_pot = c_green_pot + npp_2_greenpool
   excess_carbon = max ( 0.0_wp , c_green_pot - c_green_max )
   npp_2_greenpool = npp_2_greenpool - excess_carbon
   c_green_pot = c_green_pot - excess_carbon
   IF ( c_reserve_pot < c_reserve_optimal ) THEN
    c_reserve_pot = c_reserve_pot + excess_carbon
    npp_2_reservepool = npp_2_reservepool + excess_carbon
    excess_carbon = max ( 0.0_wp , c_reserve_pot - c_reserve_optimal )
    c_reserve_pot = c_reserve_pot - excess_carbon
    npp_2_reservepool = npp_2_reservepool - excess_carbon
   END IF
   c_litter_green_bg = c_litter_green_bg + npp_2_rootexudates
   excess_npp = excess_carbon / sec_per_day
   root_exudates = npp_2_rootexudates / sec_per_day
  END IF
  c_litter_green_ag_pot = c_litter_green_ag
  c_litter_green_bg_pot = c_litter_green_bg
  c_litter_wood_ag_pot = c_litter_wood_ag
  c_litter_wood_bg_pot = c_litter_wood_bg
  decomprate_litter_green = ( - min ( alpha_mod ** ( kappa ) * q10 ** ( (&
   topsoiltemp - referencetemp_q10 ) / 10._wp ) / tau_c_litter_green , 1._wp )&
   )
  d_c_litter_green_ag = decomprate_litter_green * c_litter_green_ag_pot
  c_litter_green_ag_pot = c_litter_green_ag_pot + d_c_litter_green_ag
  d_c_litter_green_bg = decomprate_litter_green * c_litter_green_bg_pot
  c_litter_green_bg_pot = c_litter_green_bg_pot + d_c_litter_green_bg
  d_c_litter_wood_ag = ( - min ( q10 ** ( ( topsoiltemp - referencetemp_q10 )&
   / 10._wp ) / tau_c_litter_wood , 1._wp ) * c_litter_wood_ag_pot )
  d_c_litter_wood_bg = ( - min ( q10 ** ( ( topsoiltemp - referencetemp_q10 )&
   / 10._wp ) / tau_c_litter_wood , 1._wp ) * c_litter_wood_bg_pot )
  c_litter_wood_ag_pot = c_litter_wood_ag_pot + d_c_litter_wood_ag
  c_litter_wood_bg_pot = c_litter_wood_bg_pot + d_c_litter_wood_bg
  cflx_npp_2_green_pot = npp_2_greenpool
  cflx_npp_2_wood_pot = npp_2_woodpool
  npp_2_green_via_mobilen = 0.0_wp
  npp_2_wood_via_mobilen = 0.0_wp
  cflx_npp_2_green_pot = cflx_npp_2_green_pot - npp_2_green_via_mobilen
  cflx_npp_2_wood_pot = cflx_npp_2_wood_pot - npp_2_wood_via_mobilen
  cflx_npp_2_reserve_pot = npp_2_reservepool
  cflx_littergreen_ag_2_atmos_pot = ( - fract_c_litter_green2atmos *&
   d_c_litter_green_ag )
  cflx_littergreen_bg_2_atmos_pot = ( - fract_c_litter_green2atmos *&
   d_c_litter_green_bg )
  cflx_littergreen_ag_2_slow_pot = ( - ( 1._wp - fract_c_litter_green2atmos )&
   * d_c_litter_green_ag )
  cflx_littergreen_bg_2_slow_pot = ( - ( 1._wp - fract_c_litter_green2atmos )&
   * d_c_litter_green_bg )
  cflx_litter_wood_ag_2_atmos_pot = ( - fract_c_litter_wood2atmos *&
   d_c_litter_wood_ag )
  cflx_litter_wood_bg_2_atmos_pot = ( - fract_c_litter_wood2atmos *&
   d_c_litter_wood_bg )
  cflx_litter_wood_ag_2_slow_pot = ( - ( 1._wp - fract_c_litter_wood2atmos ) *&
   d_c_litter_wood_ag )
  cflx_litter_wood_bg_2_slow_pot = ( - ( 1._wp - fract_c_litter_wood2atmos ) *&
   d_c_litter_wood_bg )
  redfact_nlimitation = 1.0_wp
  c_green = c_green + redfact_nlimitation * cflx_npp_2_green_pot
  c_woods = c_woods + redfact_nlimitation * cflx_npp_2_wood_pot
  c_reserve = c_reserve + cflx_npp_2_reserve_pot
  c_litter_green_ag = c_litter_green_ag - redfact_nlimitation * (&
   cflx_littergreen_ag_2_atmos_pot + cflx_littergreen_ag_2_slow_pot )
  c_litter_green_bg = c_litter_green_bg - redfact_nlimitation * (&
   cflx_littergreen_bg_2_atmos_pot + cflx_littergreen_bg_2_slow_pot )
  c_litter_wood_ag = c_litter_wood_ag - redfact_nlimitation * (&
   cflx_litter_wood_ag_2_atmos_pot + cflx_litter_wood_ag_2_slow_pot )
  c_litter_wood_bg = c_litter_wood_bg - redfact_nlimitation * (&
   cflx_litter_wood_bg_2_atmos_pot + cflx_litter_wood_bg_2_slow_pot )
  IF ( npp_pot_ydaymean > 0.0_wp ) THEN
   npp_act_ydaymean = ( npp_2_green_via_mobilen + npp_2_wood_via_mobilen +&
    npp_2_rootexudates + redfact_nlimitation * ( cflx_npp_2_green_pot +&
    cflx_npp_2_wood_pot ) + cflx_npp_2_reserve_pot ) / sec_per_day
  ELSE
   npp_act_ydaymean = npp_pot_ydaymean + npp_flux_correction
  END IF
  leaflitter = cflux_c_greenwood_2_litter - c_2_litter_woodpool +&
   cflx_faeces_2_lg - cflx_2_crop_harvest
  soilresp_rate = 0._wp
  CALL yasso ( yasso_io_pools ( 1 : 9 ) , weather , leaflitter , leaflit_coefv&
   , sizezero , yasso_out , fract_green_aboveground , npp_2_rootexudates ,&
   redfact_nlimitation )
  c_acid_ag1 = yasso_out ( 1 )
  c_water_ag1 = yasso_out ( 2 )
  c_ethanol_ag1 = yasso_out ( 3 )
  c_nonsoluble_ag1 = yasso_out ( 4 )
  c_acid_bg1 = yasso_out ( 5 )
  c_water_bg1 = yasso_out ( 6 )
  c_ethanol_bg1 = yasso_out ( 7 )
  c_nonsoluble_bg1 = yasso_out ( 8 )
  c_humus_1 = yasso_out ( 9 )
  soilresp_rate = yasso_out ( 10 ) / sec_per_day
  CALL yasso ( yasso_io_pools ( 10 : 18 ) , weather , c_2_litter_woodpool ,&
   woodlit_coefv , woodlittersize , yasso_out , fract_wood_aboveground ,&
   0.0_wp , redfact_nlimitation )
  c_acid_ag2 = yasso_out ( 1 )
  c_water_ag2 = yasso_out ( 2 )
  c_ethanol_ag2 = yasso_out ( 3 )
  c_nonsoluble_ag2 = yasso_out ( 4 )
  c_acid_bg2 = yasso_out ( 5 )
  c_water_bg2 = yasso_out ( 6 )
  c_ethanol_bg2 = yasso_out ( 7 )
  c_nonsoluble_bg2 = yasso_out ( 8 )
  c_humus_2 = yasso_out ( 9 )
  soilresp_rate = soilresp_rate + yasso_out ( 10 ) / sec_per_day
  cflx_crop_harvest_2_atmos = c_crop_harvest / tau_c_crop_harvest
  c_crop_harvest = c_crop_harvest - cflx_crop_harvest_2_atmos
  soilresp_rate = soilresp_rate - cflx_crop_harvest_2_atmos / sec_per_day
  cflux_herb_2_atm = ( - cflx_faeces_2_atm / sec_per_day )
  cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter / sec_per_day
  cflx_herbivory = green_2_herbivory / sec_per_day
  cflux_herb_2_littergreen = cflx_faeces_2_lg / sec_per_day
 END SUBROUTINE calc_cpools

 PURE SUBROUTINE yasso ( yasso_io_pools , weather , litter , lit_coefv ,&
  woodlittersize , yasso_out , fract_aboveground , npp_2_rootexudates ,&
  redfact_nlimit )
  USE mo_carbon_constants , ONLY: days_per_year , sec_per_year
  USE mo_jsb_physical_constants , ONLY: tmelt

  REAL ( KIND= wp ) , INTENT(IN) :: yasso_io_pools ( 1 : 9 )
  REAL ( KIND= wp ) , INTENT(IN) :: woodlittersize
  REAL ( KIND= wp ) , INTENT(IN) :: weather ( 1 : 2 )
  REAL ( KIND= wp ) , INTENT(IN) :: litter
  REAL ( KIND= wp ) , INTENT(IN) :: lit_coefv ( 1 : 5 )
  REAL ( KIND= wp ) , INTENT(IN) :: fract_aboveground
  REAL ( KIND= wp ) , INTENT(IN) :: npp_2_rootexudates
  REAL ( KIND= wp ) , INTENT(OUT) :: yasso_out ( 1 : 16 )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: redfact_nlimit
  REAL ( KIND= wp ) :: c_acid_ag
  REAL ( KIND= wp ) :: c_water_ag
  REAL ( KIND= wp ) :: c_ethanol_ag
  REAL ( KIND= wp ) :: c_nonsoluble_ag
  REAL ( KIND= wp ) :: c_acid_bg
  REAL ( KIND= wp ) :: c_water_bg
  REAL ( KIND= wp ) :: c_ethanol_bg
  REAL ( KIND= wp ) :: c_nonsoluble_bg
  REAL ( KIND= wp ) :: c_humus
  REAL ( KIND= wp ) :: cflx_litter_2_acid
  REAL ( KIND= wp ) :: cflx_litter_2_water
  REAL ( KIND= wp ) :: cflx_litter_2_ethanol
  REAL ( KIND= wp ) :: cflx_litter_2_nonsoluble
  REAL ( KIND= wp ) :: cflx_litter_2_humus
  REAL ( KIND= wp ) :: cflx_from_acid_ag
  REAL ( KIND= wp ) :: cflx_from_water_ag
  REAL ( KIND= wp ) :: cflx_from_ethanol_ag
  REAL ( KIND= wp ) :: cflx_from_nonsoluble_ag
  REAL ( KIND= wp ) :: cflx_2_acid_ag
  REAL ( KIND= wp ) :: cflx_2_water_ag
  REAL ( KIND= wp ) :: cflx_2_ethanol_ag
  REAL ( KIND= wp ) :: cflx_2_nonsoluble_ag
  REAL ( KIND= wp ) :: cflx_from_acid_bg
  REAL ( KIND= wp ) :: cflx_from_water_bg
  REAL ( KIND= wp ) :: cflx_from_ethanol_bg
  REAL ( KIND= wp ) :: cflx_from_nonsoluble_bg
  REAL ( KIND= wp ) :: cflx_from_humus
  REAL ( KIND= wp ) :: cflx_2_acid_bg
  REAL ( KIND= wp ) :: cflx_2_water_bg
  REAL ( KIND= wp ) :: cflx_2_ethanol_bg
  REAL ( KIND= wp ) :: cflx_2_nonsoluble_bg
  REAL ( KIND= wp ) :: cflx_2_humus
  REAL ( KIND= wp ) :: cflx_2_humusag
  REAL ( KIND= wp ) :: cflx_2_humusbg
  REAL ( KIND= wp ) :: soilresp_rateyasso
  REAL ( KIND= wp ) :: soilresp_ratelitterag
  REAL ( KIND= wp ) :: soilresp_ratelitterbg
  REAL ( KIND= wp ) :: d_acid
  REAL ( KIND= wp ) :: d_water
  REAL ( KIND= wp ) :: d_ethanol
  REAL ( KIND= wp ) :: d_nonsoluble
  REAL ( KIND= wp ) :: d_humus
  REAL ( KIND= wp ) :: d_litter_green
  REAL ( KIND= wp ) :: pseudo_litter_green
  REAL ( KIND= wp ) :: fract_litter_2_acid
  REAL ( KIND= wp ) :: fract_litter_2_water
  REAL ( KIND= wp ) :: fract_litter_2_ethanol
  REAL ( KIND= wp ) :: fract_litter_2_nonsoluble
  REAL ( KIND= wp ) :: fract_litter_2_humus
  REAL ( KIND= wp ) , PARAMETER :: ref_decomp_rate_acid = - 0.72_wp
  REAL ( KIND= wp ) , PARAMETER :: ref_decomp_rate_water = - 5.9_wp
  REAL ( KIND= wp ) , PARAMETER :: ref_decomp_rate_ethanol = - 0.28_wp
  REAL ( KIND= wp ) , PARAMETER :: ref_decomp_rate_nonsoluble = - 0.031_wp
  REAL ( KIND= wp ) , PARAMETER :: ref_decomp_rate_humus = - 0.0016_wp
  REAL ( KIND= wp ) , PARAMETER :: temp_p1 = 0.095_wp
  REAL ( KIND= wp ) , PARAMETER :: temp_p2 = - 0.0014_wp
  REAL ( KIND= wp ) , PARAMETER :: precip_p1 = - 1.21_wp
  REAL ( KIND= wp ) , PARAMETER :: size_p1 = - 1.71_wp
  REAL ( KIND= wp ) , PARAMETER :: size_p2 = 0.86_wp
  REAL ( KIND= wp ) , PARAMETER :: size_p3 = - 0.306_wp
  REAL ( KIND= wp ) , PARAMETER :: a_2_w = 0.99_wp
  REAL ( KIND= wp ) , PARAMETER :: a_2_e = 0.00_wp
  REAL ( KIND= wp ) , PARAMETER :: a_2_n = 0.00_wp
  REAL ( KIND= wp ) , PARAMETER :: w_2_a = 0.48_wp
  REAL ( KIND= wp ) , PARAMETER :: w_2_e = 0.00_wp
  REAL ( KIND= wp ) , PARAMETER :: w_2_n = 0.015_wp
  REAL ( KIND= wp ) , PARAMETER :: e_2_a = 0.01_wp
  REAL ( KIND= wp ) , PARAMETER :: e_2_w = 0.00_wp
  REAL ( KIND= wp ) , PARAMETER :: e_2_n = 0.95_wp
  REAL ( KIND= wp ) , PARAMETER :: n_2_a = 0.83_wp
  REAL ( KIND= wp ) , PARAMETER :: n_2_w = 0.01_wp
  REAL ( KIND= wp ) , PARAMETER :: n_2_e = 0.02_wp
  REAL ( KIND= wp ) , PARAMETER :: awen_2_h = 0.0045_wp
  REAL ( KIND= wp ) :: precip
  REAL ( KIND= wp ) :: temp
  REAL ( KIND= wp ) :: d_temp
  REAL ( KIND= wp ) :: d_precip
  REAL ( KIND= wp ) :: d_size
  REAL ( KIND= wp ) :: redfactor
  REAL ( KIND= wp ) :: dt

  dt = 1.0_wp / days_per_year
  IF ( present ( redfact_nlimit ) ) THEN
   redfactor = redfact_nlimit
  ELSE
   redfactor = 1.0_wp
  END IF
  fract_litter_2_acid = lit_coefv ( 1 )
  fract_litter_2_water = lit_coefv ( 2 )
  fract_litter_2_ethanol = lit_coefv ( 3 )
  fract_litter_2_nonsoluble = lit_coefv ( 4 )
  fract_litter_2_humus = lit_coefv ( 5 )
  c_acid_ag = yasso_io_pools ( 1 )
  c_water_ag = yasso_io_pools ( 2 )
  c_ethanol_ag = yasso_io_pools ( 3 )
  c_nonsoluble_ag = yasso_io_pools ( 4 )
  c_acid_bg = yasso_io_pools ( 5 )
  c_water_bg = yasso_io_pools ( 6 )
  c_ethanol_bg = yasso_io_pools ( 7 )
  c_nonsoluble_bg = yasso_io_pools ( 8 )
  c_humus = yasso_io_pools ( 9 )
  pseudo_litter_green = c_acid_ag + c_water_ag + c_ethanol_ag +&
   c_nonsoluble_ag + c_acid_bg + c_water_bg + c_ethanol_bg + c_nonsoluble_bg
  precip = weather ( 2 ) * sec_per_year / 1000._wp
  temp = weather ( 1 ) - tmelt
  cflx_litter_2_acid = fract_litter_2_acid * litter * days_per_year
  cflx_litter_2_water = fract_litter_2_water * litter * days_per_year
  cflx_litter_2_ethanol = fract_litter_2_ethanol * litter * days_per_year
  cflx_litter_2_nonsoluble = fract_litter_2_nonsoluble * litter *&
   days_per_year
  cflx_litter_2_humus = fract_litter_2_humus * litter * days_per_year
  d_temp = exp ( temp_p1 * temp + temp_p2 * temp ** ( 2.0_wp ) )
  d_precip = 1.0_wp - exp ( precip_p1 * precip )
  d_size = min ( 1.0_wp , ( 1.0_wp + size_p1 * woodlittersize + size_p2 *&
   woodlittersize ** ( 2.0_wp ) ) ** ( size_p3 ) )
  d_acid = redfactor * ref_decomp_rate_acid * d_temp * d_precip * d_size
  d_water = redfactor * ref_decomp_rate_water * d_temp * d_precip * d_size
  d_ethanol = redfactor * ref_decomp_rate_ethanol * d_temp * d_precip * d_size
  d_nonsoluble = redfactor * ref_decomp_rate_nonsoluble * d_temp * d_precip *&
   d_size
  d_humus = ref_decomp_rate_humus * d_temp * d_precip
  cflx_from_acid_ag = c_acid_ag * d_acid
  cflx_from_water_ag = c_water_ag * d_water
  cflx_from_ethanol_ag = c_ethanol_ag * d_ethanol
  cflx_from_nonsoluble_ag = c_nonsoluble_ag * d_nonsoluble
  cflx_from_acid_bg = c_acid_bg * d_acid
  cflx_from_water_bg = c_water_bg * d_water
  cflx_from_ethanol_bg = c_ethanol_bg * d_ethanol
  cflx_from_nonsoluble_bg = c_nonsoluble_bg * d_nonsoluble
  cflx_from_humus = c_humus * d_humus
  cflx_2_acid_ag = abs ( cflx_from_water_ag * w_2_a + cflx_from_ethanol_ag *&
   e_2_a + cflx_from_nonsoluble_ag * n_2_a )
  cflx_2_water_ag = abs ( cflx_from_acid_ag * a_2_w + cflx_from_ethanol_ag *&
   e_2_w + cflx_from_nonsoluble_ag * n_2_w )
  cflx_2_ethanol_ag = abs ( cflx_from_acid_ag * a_2_e + cflx_from_water_ag *&
   w_2_e + cflx_from_nonsoluble_ag * n_2_e )
  cflx_2_nonsoluble_ag = abs ( cflx_from_acid_ag * a_2_n + cflx_from_water_ag&
   * w_2_n + cflx_from_ethanol_ag * e_2_n )
  cflx_2_acid_bg = abs ( cflx_from_water_bg * w_2_a + cflx_from_ethanol_bg *&
   e_2_a + cflx_from_nonsoluble_bg * n_2_a )
  cflx_2_water_bg = abs ( cflx_from_acid_bg * a_2_w + cflx_from_ethanol_bg *&
   e_2_w + cflx_from_nonsoluble_bg * n_2_w )
  cflx_2_ethanol_bg = abs ( cflx_from_acid_bg * a_2_e + cflx_from_water_bg *&
   w_2_e + cflx_from_nonsoluble_bg * n_2_e )
  cflx_2_nonsoluble_bg = abs ( cflx_from_acid_bg * a_2_n + cflx_from_water_bg&
   * w_2_n + cflx_from_ethanol_bg * e_2_n )
  cflx_2_humus = abs ( cflx_from_acid_ag + cflx_from_water_ag +&
   cflx_from_ethanol_ag + cflx_from_nonsoluble_ag + cflx_from_acid_bg +&
   cflx_from_water_bg + cflx_from_ethanol_bg + cflx_from_nonsoluble_bg ) *&
   awen_2_h
  cflx_2_humusag = abs ( cflx_from_acid_ag + cflx_from_water_ag +&
   cflx_from_ethanol_ag + cflx_from_nonsoluble_ag ) * awen_2_h
  cflx_2_humusbg = abs ( cflx_from_acid_bg + cflx_from_water_bg +&
   cflx_from_ethanol_bg + cflx_from_nonsoluble_bg ) * awen_2_h
  soilresp_rateyasso = ( cflx_from_acid_ag + cflx_from_acid_bg ) * ( 1.0_wp -&
   a_2_w - a_2_e - a_2_n - awen_2_h ) + ( cflx_from_water_ag +&
   cflx_from_water_bg ) * ( 1.0_wp - w_2_a - w_2_e - w_2_n - awen_2_h ) + (&
   cflx_from_ethanol_ag + cflx_from_ethanol_bg ) * ( 1.0_wp - e_2_a - e_2_w -&
   e_2_n - awen_2_h ) + ( cflx_from_nonsoluble_ag + cflx_from_nonsoluble_bg )&
   * ( 1.0_wp - n_2_a - n_2_w - n_2_e - awen_2_h ) + cflx_from_humus
  soilresp_ratelitterag = cflx_from_acid_ag * ( 1.0_wp - a_2_w - a_2_e - a_2_n&
   - awen_2_h ) + cflx_from_water_ag * ( 1.0_wp - w_2_a - w_2_e - w_2_n -&
   awen_2_h ) + cflx_from_ethanol_ag * ( 1.0_wp - e_2_a - e_2_w - e_2_n -&
   awen_2_h ) + cflx_from_nonsoluble_ag * ( 1.0_wp - n_2_a - n_2_w - n_2_e -&
   awen_2_h )
  soilresp_ratelitterbg = cflx_from_acid_bg * ( 1.0_wp - a_2_w - a_2_e - a_2_n&
   - awen_2_h ) + cflx_from_water_bg * ( 1.0_wp - w_2_a - w_2_e - w_2_n -&
   awen_2_h ) + cflx_from_ethanol_bg * ( 1.0_wp - e_2_a - e_2_w - e_2_n -&
   awen_2_h ) + cflx_from_nonsoluble_bg * ( 1.0_wp - n_2_a - n_2_w - n_2_e -&
   awen_2_h )
  c_acid_ag = max ( 0.0_wp , c_acid_ag + ( cflx_from_acid_ag + cflx_2_acid_ag&
   + fract_aboveground * cflx_litter_2_acid ) * dt )
  c_water_ag = max ( 0.0_wp , c_water_ag + ( cflx_from_water_ag +&
   cflx_2_water_ag + fract_aboveground * cflx_litter_2_water ) * dt )
  c_ethanol_ag = max ( 0.0_wp , c_ethanol_ag + ( cflx_from_ethanol_ag +&
   cflx_2_ethanol_ag + fract_aboveground * cflx_litter_2_ethanol ) * dt )
  c_nonsoluble_ag = max ( 0.0_wp , c_nonsoluble_ag + ( cflx_from_nonsoluble_ag&
   + cflx_2_nonsoluble_ag + fract_aboveground * cflx_litter_2_nonsoluble ) *&
   dt )
  c_acid_bg = max ( 0.0_wp , c_acid_bg + ( cflx_from_acid_bg + cflx_2_acid_bg&
   + ( 1._wp - fract_aboveground ) * cflx_litter_2_acid ) * dt )
  c_water_bg = max ( 0.0_wp , c_water_bg + ( cflx_from_water_bg +&
   cflx_2_water_bg + ( 1._wp - fract_aboveground ) * cflx_litter_2_water ) *&
   dt + npp_2_rootexudates )
  c_ethanol_bg = max ( 0.0_wp , c_ethanol_bg + ( cflx_from_ethanol_bg +&
   cflx_2_ethanol_bg + ( 1._wp - fract_aboveground ) * cflx_litter_2_ethanol )&
   * dt )
  c_nonsoluble_bg = max ( 0.0_wp , c_nonsoluble_bg + ( cflx_from_nonsoluble_bg&
   + cflx_2_nonsoluble_bg + ( 1._wp - fract_aboveground ) *&
   cflx_litter_2_nonsoluble ) * dt )
  c_humus = max ( 0.0_wp , c_humus + ( cflx_from_humus + cflx_2_humus +&
   cflx_litter_2_humus ) * dt )
  IF ( pseudo_litter_green .gt. 0.0_wp ) THEN
   d_litter_green = ( c_acid_ag + c_water_ag + c_ethanol_ag + c_nonsoluble_ag&
    + c_acid_bg + c_water_bg + c_ethanol_bg + c_nonsoluble_bg ) /&
    pseudo_litter_green - 1._wp
  ELSE
   d_litter_green = 0.0_wp
  END IF
  yasso_out ( 1 ) = c_acid_ag
  yasso_out ( 2 ) = c_water_ag
  yasso_out ( 3 ) = c_ethanol_ag
  yasso_out ( 4 ) = c_nonsoluble_ag
  yasso_out ( 5 ) = c_acid_bg
  yasso_out ( 6 ) = c_water_bg
  yasso_out ( 7 ) = c_ethanol_bg
  yasso_out ( 8 ) = c_nonsoluble_bg
  yasso_out ( 9 ) = c_humus
  yasso_out ( 10 ) = soilresp_rateyasso * dt
  IF ( ( .NOT. present ( redfact_nlimit ) ) ) THEN
   yasso_out ( 11 ) = cflx_2_humusag * dt
   yasso_out ( 12 ) = cflx_2_humusbg * dt
   yasso_out ( 13 ) = ( - cflx_from_humus * dt )
   yasso_out ( 14 ) = ( - soilresp_ratelitterag * dt )
   yasso_out ( 15 ) = ( - soilresp_ratelitterbg * dt )
   yasso_out ( 16 ) = d_litter_green
  END IF
 END SUBROUTINE yasso

 SUBROUTINE relocate_carbon_fire ( burned_fract , fire_fract_wood_2_atmos ,&
  leaflit_coef , woodlit_coef , c_green , c_woods , c_reserve ,&
  cflux_dist_greenreserve_2_soil , cflux_dist_woods_2_soil , c_acid_ag1 ,&
  c_water_ag1 , c_ethanol_ag1 , c_nonsoluble_ag1 , c_acid_ag2 , c_water_ag2 ,&
  c_ethanol_ag2 , c_nonsoluble_ag2 , c_acid_bg1 , c_water_bg1 , c_ethanol_bg1&
  , c_nonsoluble_bg1 , c_acid_bg2 , c_water_bg2 , c_ethanol_bg2 ,&
  c_nonsoluble_bg2 , c_humus_1 , c_humus_2 , cflux_fire_all_2_atm )
  USE mo_carbon_constants , ONLY: fract_green_aboveground ,&
   fract_wood_aboveground

  REAL ( KIND= wp ) , INTENT(IN) :: burned_fract ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: fire_fract_wood_2_atmos
  REAL ( KIND= wp ) , INTENT(IN) :: leaflit_coef ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: woodlit_coef ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_green ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_woods ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_reserve ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: cflux_fire_all_2_atm ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: cflux_dist_greenreserve_2_soil ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: cflux_dist_woods_2_soil ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_acid_ag1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_water_ag1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_ethanol_ag1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_nonsoluble_ag1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_acid_bg1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_water_bg1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_ethanol_bg1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_nonsoluble_bg1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_humus_1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_acid_ag2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_water_ag2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_ethanol_ag2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_nonsoluble_ag2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_acid_bg2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_water_bg2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_ethanol_bg2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_nonsoluble_bg2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_humus_2 ( : )
  REAL ( KIND= wp ) :: factor ( 1 : size ( burned_fract ) )

  cflux_dist_greenreserve_2_soil ( : ) = cflux_dist_greenreserve_2_soil ( : )&
   + ( c_green ( : ) + c_reserve ( : ) ) * ( 1._wp - fract_green_aboveground )&
   * burned_fract ( : )
  cflux_dist_woods_2_soil ( : ) = cflux_dist_woods_2_soil ( : ) + c_woods ( :&
   ) * ( 1._wp - fire_fract_wood_2_atmos * fract_wood_aboveground ) *&
   burned_fract ( : )
  cflux_fire_all_2_atm ( : ) = ( ( c_green ( : ) + c_reserve ( : ) ) *&
   fract_green_aboveground + c_water_ag1 ( : ) + c_acid_ag1 ( : ) +&
   c_ethanol_ag1 ( : ) + c_nonsoluble_ag1 ( : ) + c_water_ag2 ( : ) +&
   c_acid_ag2 ( : ) + c_ethanol_ag2 ( : ) + c_nonsoluble_ag2 ( : ) + c_woods (&
   : ) * fire_fract_wood_2_atmos * fract_wood_aboveground ) * burned_fract ( :&
   )
  factor ( : ) = 1._wp - burned_fract ( : )
  c_acid_ag1 ( : ) = factor ( : ) * c_acid_ag1 ( : )
  c_water_ag1 ( : ) = factor ( : ) * c_water_ag1 ( : )
  c_ethanol_ag1 ( : ) = factor ( : ) * c_ethanol_ag1 ( : )
  c_nonsoluble_ag1 ( : ) = factor ( : ) * c_nonsoluble_ag1 ( : )
  c_acid_ag2 ( : ) = factor ( : ) * c_acid_ag2 ( : )
  c_water_ag2 ( : ) = factor ( : ) * c_water_ag2 ( : )
  c_ethanol_ag2 ( : ) = factor ( : ) * c_ethanol_ag2 ( : )
  c_nonsoluble_ag2 ( : ) = factor ( : ) * c_nonsoluble_ag2 ( : )
  factor ( : ) = ( 1._wp - fract_green_aboveground ) * burned_fract ( : ) * (&
   c_green ( : ) + c_reserve ( : ) )
  c_acid_bg1 ( : ) = c_acid_bg1 ( : ) + factor ( : ) * leaflit_coef ( 1 )
  c_water_bg1 ( : ) = c_water_bg1 ( : ) + factor ( : ) * leaflit_coef ( 2 )
  c_ethanol_bg1 ( : ) = c_ethanol_bg1 ( : ) + factor ( : ) * leaflit_coef ( 3&
   )
  c_nonsoluble_bg1 ( : ) = c_nonsoluble_bg1 ( : ) + factor ( : ) *&
   leaflit_coef ( 4 )
  c_humus_1 ( : ) = c_humus_1 ( : ) + factor ( : ) * leaflit_coef ( 5 )
  factor ( : ) = ( 1._wp - fire_fract_wood_2_atmos ) * fract_wood_aboveground&
   * burned_fract ( : ) * c_woods ( : )
  c_acid_ag2 ( : ) = c_acid_ag2 ( : ) + factor ( : ) * woodlit_coef ( 1 )
  c_water_ag2 ( : ) = c_water_ag2 ( : ) + factor ( : ) * woodlit_coef ( 2 )
  c_ethanol_ag2 ( : ) = c_ethanol_ag2 ( : ) + factor ( : ) * woodlit_coef ( 3&
   )
  c_nonsoluble_ag2 ( : ) = c_nonsoluble_ag2 ( : ) + factor ( : ) *&
   woodlit_coef ( 4 )
  c_humus_2 ( : ) = c_humus_2 ( : ) + factor ( : ) * woodlit_coef ( 5 )
  factor ( : ) = ( 1._wp - fract_wood_aboveground ) * burned_fract ( : ) *&
   c_woods ( : )
  c_acid_bg2 ( : ) = c_acid_bg2 ( : ) + factor ( : ) * woodlit_coef ( 1 )
  c_water_bg2 ( : ) = c_water_bg2 ( : ) + factor ( : ) * woodlit_coef ( 2 )
  c_ethanol_bg2 ( : ) = c_ethanol_bg2 ( : ) + factor ( : ) * woodlit_coef ( 3&
   )
  c_nonsoluble_bg2 ( : ) = c_nonsoluble_bg2 ( : ) + factor ( : ) *&
   woodlit_coef ( 4 )
  c_humus_2 ( : ) = c_humus_2 ( : ) + factor ( : ) * woodlit_coef ( 5 )
  factor ( : ) = 1._wp - burned_fract ( : )
  c_green ( : ) = c_green ( : ) * factor ( : )
  c_woods ( : ) = c_woods ( : ) * factor ( : )
  c_reserve ( : ) = c_reserve ( : ) * factor ( : )
 END SUBROUTINE relocate_carbon_fire

 SUBROUTINE relocate_carbon_damage ( damaged_fract , leaflit_coef ,&
  woodlit_coef , c_green , c_woods , c_reserve ,&
  cflux_dist_greenreserve_2_soil , cflux_dist_woods_2_soil , c_acid_ag1 ,&
  c_water_ag1 , c_ethanol_ag1 , c_nonsoluble_ag1 , c_acid_ag2 , c_water_ag2 ,&
  c_ethanol_ag2 , c_nonsoluble_ag2 , c_acid_bg1 , c_water_bg1 , c_ethanol_bg1&
  , c_nonsoluble_bg1 , c_acid_bg2 , c_water_bg2 , c_ethanol_bg2 ,&
  c_nonsoluble_bg2 , c_humus_1 , c_humus_2 )
  USE mo_carbon_constants , ONLY: fract_green_aboveground ,&
   fract_wood_aboveground

  REAL ( KIND= wp ) , INTENT(IN) :: damaged_fract ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: leaflit_coef ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: woodlit_coef ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_green ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_woods ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_reserve ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: cflux_dist_greenreserve_2_soil ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: cflux_dist_woods_2_soil ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_acid_ag1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_water_ag1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_ethanol_ag1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_nonsoluble_ag1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_acid_bg1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_water_bg1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_ethanol_bg1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_nonsoluble_bg1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_humus_1 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_acid_ag2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_water_ag2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_ethanol_ag2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_nonsoluble_ag2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_acid_bg2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_water_bg2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_ethanol_bg2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_nonsoluble_bg2 ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: c_humus_2 ( : )
  REAL ( KIND= wp ) :: factor ( 1 : size ( damaged_fract ) )

  cflux_dist_woods_2_soil ( : ) = cflux_dist_woods_2_soil ( : ) + c_woods ( :&
   ) * damaged_fract ( : )
  cflux_dist_greenreserve_2_soil ( : ) = cflux_dist_greenreserve_2_soil ( : )&
   + ( c_green ( : ) + c_reserve ( : ) ) * damaged_fract ( : )
  factor ( : ) = ( c_green ( : ) + c_reserve ( : ) ) * fract_green_aboveground&
   * damaged_fract ( : )
  c_acid_ag1 ( : ) = c_acid_ag1 ( : ) + factor ( : ) * leaflit_coef ( 1 )
  c_water_ag1 ( : ) = c_water_ag1 ( : ) + factor ( : ) * leaflit_coef ( 2 )
  c_ethanol_ag1 ( : ) = c_ethanol_ag1 ( : ) + factor ( : ) * leaflit_coef ( 3&
   )
  c_nonsoluble_ag1 ( : ) = c_nonsoluble_ag1 ( : ) + factor ( : ) *&
   leaflit_coef ( 4 )
  factor ( : ) = ( c_green ( : ) + c_reserve ( : ) ) * ( 1._wp -&
   fract_green_aboveground ) * damaged_fract ( : )
  c_acid_bg1 ( : ) = c_acid_bg1 ( : ) + factor ( : ) * leaflit_coef ( 1 )
  c_water_bg1 ( : ) = c_water_bg1 ( : ) + factor ( : ) * leaflit_coef ( 2 )
  c_ethanol_bg1 ( : ) = c_ethanol_bg1 ( : ) + factor ( : ) * leaflit_coef ( 3&
   )
  c_nonsoluble_bg1 ( : ) = c_nonsoluble_bg1 ( : ) + factor ( : ) *&
   leaflit_coef ( 4 )
  factor ( : ) = ( c_green ( : ) + c_reserve ( : ) ) * damaged_fract ( : )
  c_humus_1 ( : ) = c_humus_1 ( : ) + factor ( : ) * leaflit_coef ( 5 )
  factor ( : ) = fract_wood_aboveground * damaged_fract ( : ) * c_woods ( : )
  c_acid_ag2 ( : ) = c_acid_ag2 ( : ) + factor ( : ) * woodlit_coef ( 1 )
  c_water_ag2 ( : ) = c_water_ag2 ( : ) + factor ( : ) * woodlit_coef ( 2 )
  c_ethanol_ag2 ( : ) = c_ethanol_ag2 ( : ) + factor ( : ) * woodlit_coef ( 3&
   )
  c_nonsoluble_ag2 ( : ) = c_nonsoluble_ag2 ( : ) + factor ( : ) *&
   woodlit_coef ( 4 )
  factor ( : ) = ( 1._wp - fract_wood_aboveground ) * damaged_fract ( : ) *&
   c_woods ( : )
  c_acid_bg2 ( : ) = c_acid_bg2 ( : ) + factor ( : ) * woodlit_coef ( 1 )
  c_water_bg2 ( : ) = c_water_bg2 ( : ) + factor ( : ) * woodlit_coef ( 2 )
  c_ethanol_bg2 ( : ) = c_ethanol_bg2 ( : ) + factor ( : ) * woodlit_coef ( 3&
   )
  c_nonsoluble_bg2 ( : ) = c_nonsoluble_bg2 ( : ) + factor ( : ) *&
   woodlit_coef ( 4 )
  factor ( : ) = c_woods ( : ) * damaged_fract ( : )
  c_humus_2 ( : ) = c_humus_2 ( : ) + factor ( : ) * woodlit_coef ( 5 )
  factor ( : ) = 1._wp - damaged_fract ( : )
  c_green ( : ) = c_green ( : ) * factor ( : )
  c_woods ( : ) = c_woods ( : ) * factor ( : )
  c_reserve ( : ) = c_reserve ( : ) * factor ( : )
 END SUBROUTINE relocate_carbon_damage

 SUBROUTINE get_per_tile ( cpool_per_tile , cpool_per_canopy ,&
  veg_fract_correction , fract_fpc_max )

  REAL ( KIND= wp ) , INTENT(IN) :: fract_fpc_max ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: veg_fract_correction ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: cpool_per_canopy ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: cpool_per_tile ( : )
  INTEGER :: jsb_hori

!$acc data &
!$acc present(fract_fpc_max,veg_fract_correction,cpool_per_canopy &
!$acc ,cpool_per_tile)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( cpool_per_canopy , 1 ) , 1
   cpool_per_tile ( jsb_hori ) = cpool_per_canopy ( jsb_hori ) *&
    veg_fract_correction ( jsb_hori ) * fract_fpc_max ( jsb_hori )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE get_per_tile

END MODULE mo_carbon_process

