MODULE mo_carbon_interface
 USE mo_jsb_control , ONLY: debug_on
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message , finish
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_grid_class , ONLY: t_jsb_grid
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract , t_jsb_aggregator
 USE mo_jsb_process_class , ONLY: t_jsb_process
 USE mo_jsb_task_class , ONLY: t_jsb_process_task , t_jsb_task_options
 USE mo_jsb_process_class , ONLY: carbon_ , assimi_ , pheno_ , disturb_ , a2l_&
  , l2a_
 USE mo_pheno_config_class , ONLY: t_pheno_config
 USE mo_a2l_memory_class , ONLY: t_a2l_memory
 USE mo_assimi_memory_class , ONLY: t_assimi_memory
 USE mo_pheno_memory_class , ONLY: t_pheno_memory
 USE mo_carbon_memory_class , ONLY: t_carbon_memory
 USE mo_l2a_memory_class , ONLY: t_l2a_memory
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_carbon_interface"
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_c_npp_pot_allocation
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_c_npp_pot_allocation
  PROCEDURE , NOPASS :: aggregate => aggregate_c_npp_pot_allocation
 END TYPE tsk_c_npp_pot_allocation
 PRIVATE :: create_task_c_npp_pot_allocation
 PRIVATE :: update_c_npp_pot_allocation
 PRIVATE :: aggregate_c_npp_pot_allocation
 PUBLIC :: register_carbon_tasks
 PUBLIC :: recalc_per_tile_vars
 PUBLIC :: calculate_current_c_ag_1_and_bg_sums
 PUBLIC :: rescale_carbon_upon_reference_area_change
 PUBLIC :: calculate_current_c_ta_state_sum
 PUBLIC :: check_carbon_conservation
 PUBLIC :: yday_carbon_conservation_test
 INTERFACE tsk_c_npp_pot_allocation
  PROCEDURE create_task_c_npp_pot_allocation
 END INTERFACE

CONTAINS
 FUNCTION create_task_c_npp_pot_allocation ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_c_npp_pot_allocation :: return_ptr )
  CALL return_ptr % construct ( name ="C_NPP_pot_allocation" , process_id =&
   carbon_ , owner_model_id = model_id )
 END FUNCTION create_task_c_npp_pot_allocation

 SUBROUTINE register_carbon_tasks ( this , model_id )

  CLASS ( t_jsb_process ) , INTENT(INOUT) :: this
  INTEGER , INTENT(IN) :: model_id

  CALL this % register_task ( tsk_c_npp_pot_allocation ( model_id ) )
 END SUBROUTINE register_carbon_tasks

 SUBROUTINE update_c_npp_pot_allocation ( tile , options )
  USE mo_carbon_process , ONLY: calc_cpools
  USE mo_jsb_time , ONLY: is_newday , is_newyear , time_steps_per_day ,&
   is_time_experiment_start
  USE mo_carbon_constants , ONLY: molarmassco2_kg , sec_per_day

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  LOGICAL :: lstart
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) , ALLOCATABLE :: precip_total ( : )
  REAL ( KIND= wp ) :: maxlai ( 1 : options % nc )
  REAL ( KIND= wp ) :: old_c_state_sum ( 1 : options % nc )
  REAL ( KIND= wp ) :: current_fluxes ( 1 : options % nc )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_carbon_interface" //&
  ":update_C_NPP_pot_allocation"
  TYPE ( t_pheno_config ) , POINTER :: pheno__conf
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  REAL ( KIND= wp ) , POINTER :: pseudo_temp ( : )
  REAL ( KIND= wp ) , POINTER :: n_pseudo_temp ( : )
  REAL ( KIND= wp ) , POINTER :: f_pseudo_temp ( : )
  REAL ( KIND= wp ) , POINTER :: pseudo_temp_yday ( : )
  REAL ( KIND= wp ) , POINTER :: pseudo_precip ( : )
  REAL ( KIND= wp ) , POINTER :: pseudo_precip_yday ( : )
  REAL ( KIND= wp ) , POINTER :: n_pseudo_precip ( : )
  REAL ( KIND= wp ) , POINTER :: f_pseudo_precip ( : )
  REAL ( KIND= wp ) , POINTER :: topsoiltemp ( : )
  REAL ( KIND= wp ) , POINTER :: alpha ( : )
  REAL ( KIND= wp ) , POINTER :: cconservation_calccpools ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_c_greenwood_2_litter ( : )
  REAL ( KIND= wp ) , POINTER :: c_green ( : )
  REAL ( KIND= wp ) , POINTER :: c_reserve ( : )
  REAL ( KIND= wp ) , POINTER :: c_woods ( : )
  REAL ( KIND= wp ) , POINTER :: c_crop_harvest ( : )
  REAL ( KIND= wp ) , POINTER :: c_litter_green_ag ( : )
  REAL ( KIND= wp ) , POINTER :: c_litter_green_bg ( : )
  REAL ( KIND= wp ) , POINTER :: c_litter_wood_ag ( : )
  REAL ( KIND= wp ) , POINTER :: c_litter_wood_bg ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_2 ( : )
  REAL ( KIND= wp ) , POINTER :: gross_assimilation_ca ( : )
  REAL ( KIND= wp ) , POINTER :: npp_pot_rate_ca ( : )
  REAL ( KIND= wp ) , POINTER :: lai_sum ( : )
  REAL ( KIND= wp ) , POINTER :: npp_pot_sum ( : )
  REAL ( KIND= wp ) , POINTER :: gpp_sum ( : )
  REAL ( KIND= wp ) , POINTER :: topsoiltemp_sum ( : )
  REAL ( KIND= wp ) , POINTER :: alpha_sum ( : )
  REAL ( KIND= wp ) , POINTER :: lai_ydaymean ( : )
  REAL ( KIND= wp ) , POINTER :: lai_yydaymean ( : )
  REAL ( KIND= wp ) , POINTER :: npp_pot_ydaymean ( : )
  REAL ( KIND= wp ) , POINTER :: npp_act_ydaymean ( : )
  REAL ( KIND= wp ) , POINTER :: gpp_ydaymean ( : )
  REAL ( KIND= wp ) , POINTER :: topsoiltemp_ydaymean ( : )
  REAL ( KIND= wp ) , POINTER :: alpha_ydaymean ( : )
  REAL ( KIND= wp ) , POINTER :: soil_respiration ( : )
  REAL ( KIND= wp ) , POINTER :: npp_flux_correction ( : )
  REAL ( KIND= wp ) , POINTER :: excess_npp ( : )
  REAL ( KIND= wp ) , POINTER :: root_exudates ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_c_green_2_herb ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_herb_2_littergreen ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_herb_2_atm ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_npp_2_atm_yday_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_npp_2_atm_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_soilresp_2_atm_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_herb_2_atm_ta ( : )
  REAL ( KIND= wp ) , POINTER :: veg_fract_correction ( : )
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc_max ( : )
  REAL ( KIND= wp ) , POINTER :: t_air ( : )
  REAL ( KIND= wp ) , POINTER :: rain ( : )
  REAL ( KIND= wp ) , POINTER :: snow ( : )

  IF ( ( .NOT. tile % is_process_active ( carbon_ ) ) ) THEN
   RETURN
  END IF
  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_carbon_interface" //":update_C_NPP_pot_allocation"&
    ) ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  SELECT TYPE ( cc => model % processes ( pheno_ ) % p % config )
  TYPE IS ( t_pheno_config )
   pheno__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
  TYPE IS ( t_carbon_memory )
   carbon__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  pseudo_temp => carbon__mem % pseudo_temp % ptr ( ics : ice , iblk )
  pseudo_temp_yday => carbon__mem % pseudo_temp_yday % ptr ( ics : ice , iblk&
   )
  n_pseudo_temp => carbon__mem % n_pseudo_temp % ptr ( ics : ice , iblk )
  f_pseudo_temp => carbon__mem % f_pseudo_temp % ptr ( ics : ice , iblk )
  pseudo_precip => carbon__mem % pseudo_precip % ptr ( ics : ice , iblk )
  pseudo_precip_yday => carbon__mem % pseudo_precip_yday % ptr ( ics : ice ,&
   iblk )
  n_pseudo_precip => carbon__mem % n_pseudo_precip % ptr ( ics : ice , iblk )
  f_pseudo_precip => carbon__mem % f_pseudo_precip % ptr ( ics : ice , iblk )
  topsoiltemp => carbon__mem % topsoiltemp % ptr ( ics : ice , iblk )
  alpha => carbon__mem % alpha % ptr ( ics : ice , iblk )
  cconservation_calccpools => carbon__mem % cconservation_calccpools % ptr (&
   ics : ice , iblk )
  cflux_c_greenwood_2_litter => carbon__mem % cflux_c_greenwood_2_litter % ptr&
   ( ics : ice , iblk )
  c_green => carbon__mem % c_green % ptr ( ics : ice , iblk )
  c_reserve => carbon__mem % c_reserve % ptr ( ics : ice , iblk )
  c_woods => carbon__mem % c_woods % ptr ( ics : ice , iblk )
  c_crop_harvest => carbon__mem % c_crop_harvest % ptr ( ics : ice , iblk )
  c_litter_green_ag => carbon__mem % c_litter_green_ag % ptr ( ics : ice ,&
   iblk )
  c_litter_green_bg => carbon__mem % c_litter_green_bg % ptr ( ics : ice ,&
   iblk )
  c_litter_wood_ag => carbon__mem % c_litter_wood_ag % ptr ( ics : ice , iblk&
   )
  c_litter_wood_bg => carbon__mem % c_litter_wood_bg % ptr ( ics : ice , iblk&
   )
  c_acid_ag1 => carbon__mem % c_acid_ag1 % ptr ( ics : ice , iblk )
  c_water_ag1 => carbon__mem % c_water_ag1 % ptr ( ics : ice , iblk )
  c_ethanol_ag1 => carbon__mem % c_ethanol_ag1 % ptr ( ics : ice , iblk )
  c_nonsoluble_ag1 => carbon__mem % c_nonsoluble_ag1 % ptr ( ics : ice , iblk&
   )
  c_acid_bg1 => carbon__mem % c_acid_bg1 % ptr ( ics : ice , iblk )
  c_water_bg1 => carbon__mem % c_water_bg1 % ptr ( ics : ice , iblk )
  c_ethanol_bg1 => carbon__mem % c_ethanol_bg1 % ptr ( ics : ice , iblk )
  c_nonsoluble_bg1 => carbon__mem % c_nonsoluble_bg1 % ptr ( ics : ice , iblk&
   )
  c_humus_1 => carbon__mem % c_humus_1 % ptr ( ics : ice , iblk )
  c_acid_ag2 => carbon__mem % c_acid_ag2 % ptr ( ics : ice , iblk )
  c_water_ag2 => carbon__mem % c_water_ag2 % ptr ( ics : ice , iblk )
  c_ethanol_ag2 => carbon__mem % c_ethanol_ag2 % ptr ( ics : ice , iblk )
  c_nonsoluble_ag2 => carbon__mem % c_nonsoluble_ag2 % ptr ( ics : ice , iblk&
   )
  c_acid_bg2 => carbon__mem % c_acid_bg2 % ptr ( ics : ice , iblk )
  c_water_bg2 => carbon__mem % c_water_bg2 % ptr ( ics : ice , iblk )
  c_ethanol_bg2 => carbon__mem % c_ethanol_bg2 % ptr ( ics : ice , iblk )
  c_nonsoluble_bg2 => carbon__mem % c_nonsoluble_bg2 % ptr ( ics : ice , iblk&
   )
  c_humus_2 => carbon__mem % c_humus_2 % ptr ( ics : ice , iblk )
  gross_assimilation_ca => assimi__mem % gross_assimilation_ca % ptr ( ics :&
   ice , iblk )
  npp_pot_rate_ca => assimi__mem % npp_pot_rate_ca % ptr ( ics : ice , iblk )
  lai_sum => carbon__mem % lai_sum % ptr ( ics : ice , iblk )
  npp_pot_sum => carbon__mem % npp_pot_sum % ptr ( ics : ice , iblk )
  gpp_sum => carbon__mem % gpp_sum % ptr ( ics : ice , iblk )
  topsoiltemp_sum => carbon__mem % topsoiltemp_sum % ptr ( ics : ice , iblk )
  alpha_sum => carbon__mem % alpha_sum % ptr ( ics : ice , iblk )
  lai_ydaymean => carbon__mem % lai_ydaymean % ptr ( ics : ice , iblk )
  lai_yydaymean => carbon__mem % lai_yydaymean % ptr ( ics : ice , iblk )
  npp_pot_ydaymean => carbon__mem % npp_pot_ydaymean % ptr ( ics : ice , iblk&
   )
  npp_act_ydaymean => carbon__mem % npp_act_ydaymean % ptr ( ics : ice , iblk&
   )
  gpp_ydaymean => carbon__mem % gpp_ydaymean % ptr ( ics : ice , iblk )
  topsoiltemp_ydaymean => carbon__mem % topsoiltemp_ydaymean % ptr ( ics : ice&
   , iblk )
  alpha_ydaymean => carbon__mem % alpha_ydaymean % ptr ( ics : ice , iblk )
  soil_respiration => carbon__mem % soil_respiration % ptr ( ics : ice , iblk&
   )
  npp_flux_correction => carbon__mem % npp_flux_correction % ptr ( ics : ice ,&
   iblk )
  excess_npp => carbon__mem % excess_npp % ptr ( ics : ice , iblk )
  root_exudates => carbon__mem % root_exudates % ptr ( ics : ice , iblk )
  cflux_c_green_2_herb => carbon__mem % cflux_c_green_2_herb % ptr ( ics : ice&
   , iblk )
  cflux_herb_2_littergreen => carbon__mem % cflux_herb_2_littergreen % ptr (&
   ics : ice , iblk )
  cflux_herb_2_atm => carbon__mem % cflux_herb_2_atm % ptr ( ics : ice , iblk&
   )
  cflux_npp_2_atm_yday_ta => carbon__mem % cflux_npp_2_atm_yday_ta % ptr ( ics&
   : ice , iblk )
  cflux_npp_2_atm_ta => carbon__mem % cflux_npp_2_atm_ta % ptr ( ics : ice ,&
   iblk )
  cflux_soilresp_2_atm_ta => carbon__mem % cflux_soilresp_2_atm_ta % ptr ( ics&
   : ice , iblk )
  cflux_herb_2_atm_ta => carbon__mem % cflux_herb_2_atm_ta % ptr ( ics : ice ,&
   iblk )
  veg_fract_correction => pheno__mem % veg_fract_correction % ptr ( ics : ice&
   , iblk )
  lai => pheno__mem % lai % ptr ( ics : ice , iblk )
  fract_fpc_max => pheno__mem % fract_fpc_max % ptr ( ics : ice , iblk )
  t_air => a2l__mem % t_air % ptr ( ics : ice , iblk )
  rain => a2l__mem % rain % ptr ( ics : ice , iblk )
  snow => a2l__mem % snow % ptr ( ics : ice , iblk )
  IF ( debug_on ( ) .AND. iblk == 1 .AND. is_newday ( ) ) THEN
   CALL message ( trim ("mo_carbon_interface" //":update_C_NPP_pot_allocation"&
    ) ,"update_C_NPP_pot_allocation: First time step of new day" )
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 .AND. is_newyear ( ) ) THEN
   CALL message ( trim ("mo_carbon_interface" //":update_C_NPP_pot_allocation"&
    ) ,"update_C_NPP_pot_allocation: First time step of new year" )
  END IF
  cflux_npp_2_atm_yday_ta = 0._wp
  cflux_npp_2_atm_ta = 0._wp
  cflux_soilresp_2_atm_ta = 0._wp
  cflux_herb_2_atm_ta = 0._wp
  pseudo_temp = t_air * n_pseudo_temp + f_pseudo_temp * pseudo_temp
  precip_total = rain + snow
  pseudo_precip = precip_total * n_pseudo_precip + f_pseudo_precip *&
   pseudo_precip
  lstart = is_time_experiment_start ( )
  IF ( ( .NOT. is_newday ( ) ) .OR. lstart ) THEN
   lai_sum = lai_sum + lai
   npp_pot_sum = npp_pot_sum + npp_pot_rate_ca
   gpp_sum = gpp_sum + gross_assimilation_ca
   topsoiltemp_sum = topsoiltemp_sum + topsoiltemp
   alpha_sum = alpha_sum + alpha
  ELSE
   lai_yydaymean = lai_ydaymean
   lai_ydaymean = lai_sum / time_steps_per_day ( )
   npp_pot_ydaymean = npp_pot_sum / time_steps_per_day ( )
   gpp_ydaymean = gpp_sum / time_steps_per_day ( )
   topsoiltemp_ydaymean = topsoiltemp_sum / time_steps_per_day ( )
   alpha_ydaymean = alpha_sum / time_steps_per_day ( )
   lai_sum = lai
   npp_pot_sum = npp_pot_rate_ca
   gpp_sum = gross_assimilation_ca
   topsoiltemp_sum = topsoiltemp
   alpha_sum = alpha
   pseudo_temp_yday = pseudo_temp
   pseudo_precip_yday = pseudo_precip
   IF ( pheno__conf % l_forestregrowth ) THEN
    maxlai ( : ) = pheno__mem % maxlai_allom % ptr ( ics : ice , iblk )
   ELSE
    maxlai ( : ) = model % lctlib ( tile % lcts ( 1 ) % lib_id ) % maxlai
   END IF
   CALL calculate_current_c_ta_state_sum ( tile , options , old_c_state_sum (&
    : ) )
   CALL calc_cpools ( lai_ydaymean , lai_yydaymean , maxlai , npp_pot_ydaymean&
    , topsoiltemp_ydaymean , alpha_ydaymean , model % lctlib ( tile % lcts ( 1&
    ) % lib_id ) % fract_npp_2_woodpool , model % lctlib ( tile % lcts ( 1 ) %&
    lib_id ) % fract_npp_2_reservepool , model % lctlib ( tile % lcts ( 1 ) %&
    lib_id ) % fract_npp_2_exudates , model % lctlib ( tile % lcts ( 1 ) %&
    lib_id ) % fract_green_2_herbivory , model % lctlib ( tile % lcts ( 1 ) %&
    lib_id ) % tau_c_litter_leaf , model % lctlib ( tile % lcts ( 1 ) % lib_id&
    ) % tau_c_litter_wood , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    tau_c_woods , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    lai_shed_constant , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    fract_c_litter_green2atmos , model % lctlib ( tile % lcts ( 1 ) % lib_id )&
    % max_c_content_woods , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    specificleafarea_c , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    reservec2leafc , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    phenologytype , c_green , c_woods , c_reserve , c_crop_harvest ,&
    c_litter_green_ag , c_litter_green_bg , c_litter_wood_ag ,&
    c_litter_wood_bg , soil_respiration , npp_flux_correction , excess_npp ,&
    root_exudates , cflux_c_greenwood_2_litter , cflux_c_green_2_herb ,&
    cflux_herb_2_littergreen , cflux_herb_2_atm , npp_act_ydaymean , temp2_30d&
    = pseudo_temp_yday , precip_30d = pseudo_precip_yday , c_acid_ag1 =&
    c_acid_ag1 , c_water_ag1 = c_water_ag1 , c_ethanol_ag1 = c_ethanol_ag1 ,&
    c_nonsoluble_ag1 = c_nonsoluble_ag1 , c_acid_bg1 = c_acid_bg1 ,&
    c_water_bg1 = c_water_bg1 , c_ethanol_bg1 = c_ethanol_bg1 ,&
    c_nonsoluble_bg1 = c_nonsoluble_bg1 , c_humus_1 = c_humus_1 , c_acid_ag2 =&
    c_acid_ag2 , c_water_ag2 = c_water_ag2 , c_ethanol_ag2 = c_ethanol_ag2 ,&
    c_nonsoluble_ag2 = c_nonsoluble_ag2 , c_acid_bg2 = c_acid_bg2 ,&
    c_water_bg2 = c_water_bg2 , c_ethanol_bg2 = c_ethanol_bg2 ,&
    c_nonsoluble_bg2 = c_nonsoluble_bg2 , c_humus_2 = c_humus_2 ,&
    leaflit_coef1 = model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    leaflit_coef ( 1 ) , leaflit_coef2 = model % lctlib ( tile % lcts ( 1 ) %&
    lib_id ) % leaflit_coef ( 2 ) , leaflit_coef3 = model % lctlib ( tile %&
    lcts ( 1 ) % lib_id ) % leaflit_coef ( 3 ) , leaflit_coef4 = model %&
    lctlib ( tile % lcts ( 1 ) % lib_id ) % leaflit_coef ( 4 ) , leaflit_coef5&
    = model % lctlib ( tile % lcts ( 1 ) % lib_id ) % leaflit_coef ( 5 ) ,&
    woodlit_coef1 = model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    woodlit_coef ( 1 ) , woodlit_coef2 = model % lctlib ( tile % lcts ( 1 ) %&
    lib_id ) % woodlit_coef ( 2 ) , woodlit_coef3 = model % lctlib ( tile %&
    lcts ( 1 ) % lib_id ) % woodlit_coef ( 3 ) , woodlit_coef4 = model %&
    lctlib ( tile % lcts ( 1 ) % lib_id ) % woodlit_coef ( 4 ) , woodlit_coef5&
    = model % lctlib ( tile % lcts ( 1 ) % lib_id ) % woodlit_coef ( 5 ) ,&
    woodlittersize = model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    woodlittersize )
   CALL calculate_current_c_ag_1_and_bg_sums ( tile , options )
   CALL recalc_per_tile_vars ( tile , options , c_green = c_green ( : ) ,&
    c_woods = c_woods ( : ) , c_reserve = c_reserve ( : ) , c_crop_harvest =&
    c_crop_harvest ( : ) , c_acid_ag1 = c_acid_ag1 ( : ) , c_water_ag1 =&
    c_water_ag1 ( : ) , c_ethanol_ag1 = c_ethanol_ag1 ( : ) , c_nonsoluble_ag1&
    = c_nonsoluble_ag1 ( : ) , c_acid_ag2 = c_acid_ag2 ( : ) , c_water_ag2 =&
    c_water_ag2 ( : ) , c_ethanol_ag2 = c_ethanol_ag2 ( : ) , c_nonsoluble_ag2&
    = c_nonsoluble_ag2 ( : ) , c_acid_bg1 = c_acid_bg1 ( : ) , c_water_bg1 =&
    c_water_bg1 ( : ) , c_ethanol_bg1 = c_ethanol_bg1 ( : ) , c_nonsoluble_bg1&
    = c_nonsoluble_bg1 ( : ) , c_acid_bg2 = c_acid_bg2 ( : ) , c_water_bg2 =&
    c_water_bg2 ( : ) , c_ethanol_bg2 = c_ethanol_bg2 ( : ) , c_nonsoluble_bg2&
    = c_nonsoluble_bg2 ( : ) , c_humus_1 = c_humus_1 ( : ) , c_humus_2 =&
    c_humus_2 ( : ) , root_exudates = root_exudates ( : ) , soil_respiration =&
    soil_respiration ( : ) , npp_flux_correction = npp_flux_correction ( : ) ,&
    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter ( : ) ,&
    cflux_c_green_2_herb = cflux_c_green_2_herb ( : ) , cflux_herb_2_atm =&
    cflux_herb_2_atm ( : ) , npp_act_ydaymean = npp_act_ydaymean ( : ) ,&
    npp_pot_ydaymean = npp_pot_ydaymean ( : ) , gpp_ydaymean = gpp_ydaymean (&
    : ) )
   current_fluxes = ( npp_act_ydaymean + cflux_herb_2_atm + soil_respiration )&
    * sec_per_day * veg_fract_correction * fract_fpc_max
   CALL check_carbon_conservation ( tile , options , old_c_state_sum ( : ) ,&
    current_fluxes ( : ) , cconservation_calccpools ( : ) )
  END IF
  cflux_npp_2_atm_ta = molarmassco2_kg * ( ( - veg_fract_correction *&
   fract_fpc_max * ( npp_pot_rate_ca - ( npp_pot_ydaymean - npp_act_ydaymean )&
   ) ) )
  cflux_soilresp_2_atm_ta = molarmassco2_kg * ( ( - veg_fract_correction *&
   fract_fpc_max * soil_respiration ) )
  cflux_herb_2_atm_ta = molarmassco2_kg * ( ( - veg_fract_correction *&
   fract_fpc_max * cflux_herb_2_atm ) )
  cflux_npp_2_atm_yday_ta = molarmassco2_kg * ( ( - veg_fract_correction *&
   fract_fpc_max * npp_act_ydaymean ) )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_carbon_interface" //":update_C_NPP_pot_allocation"&
    ) ,"Finished." )
  END IF
 END SUBROUTINE update_c_npp_pot_allocation

 SUBROUTINE calculate_current_c_ta_state_sum ( tile , options ,&
  current_c_state_sum )

  CLASS ( t_jsb_tile_abstract ) , TARGET , INTENT(IN) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  REAL ( KIND= wp ) , INTENT(INOUT) :: current_c_state_sum ( : )
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_carbon_interface" //&
  ":calculate_current_c_ta_state_sum"
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: iblk
  REAL ( KIND= wp ) , POINTER :: c_green_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_reserve_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_woods_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_crop_harvest_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_2_ta ( : )
  INTEGER :: claw_induction_2

  SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
  TYPE IS ( t_carbon_memory )
   carbon__mem => mm
  END SELECT
  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  c_green_ta => carbon__mem % c_green_ta % ptr ( ics : ice , iblk )
  c_reserve_ta => carbon__mem % c_reserve_ta % ptr ( ics : ice , iblk )
  c_woods_ta => carbon__mem % c_woods_ta % ptr ( ics : ice , iblk )
  c_crop_harvest_ta => carbon__mem % c_crop_harvest_ta % ptr ( ics : ice ,&
   iblk )
  c_acid_ag1_ta => carbon__mem % c_acid_ag1_ta % ptr ( ics : ice , iblk )
  c_water_ag1_ta => carbon__mem % c_water_ag1_ta % ptr ( ics : ice , iblk )
  c_ethanol_ag1_ta => carbon__mem % c_ethanol_ag1_ta % ptr ( ics : ice , iblk&
   )
  c_nonsoluble_ag1_ta => carbon__mem % c_nonsoluble_ag1_ta % ptr ( ics : ice ,&
   iblk )
  c_acid_bg1_ta => carbon__mem % c_acid_bg1_ta % ptr ( ics : ice , iblk )
  c_water_bg1_ta => carbon__mem % c_water_bg1_ta % ptr ( ics : ice , iblk )
  c_ethanol_bg1_ta => carbon__mem % c_ethanol_bg1_ta % ptr ( ics : ice , iblk&
   )
  c_nonsoluble_bg1_ta => carbon__mem % c_nonsoluble_bg1_ta % ptr ( ics : ice ,&
   iblk )
  c_humus_1_ta => carbon__mem % c_humus_1_ta % ptr ( ics : ice , iblk )
  c_acid_ag2_ta => carbon__mem % c_acid_ag2_ta % ptr ( ics : ice , iblk )
  c_water_ag2_ta => carbon__mem % c_water_ag2_ta % ptr ( ics : ice , iblk )
  c_ethanol_ag2_ta => carbon__mem % c_ethanol_ag2_ta % ptr ( ics : ice , iblk&
   )
  c_nonsoluble_ag2_ta => carbon__mem % c_nonsoluble_ag2_ta % ptr ( ics : ice ,&
   iblk )
  c_acid_bg2_ta => carbon__mem % c_acid_bg2_ta % ptr ( ics : ice , iblk )
  c_water_bg2_ta => carbon__mem % c_water_bg2_ta % ptr ( ics : ice , iblk )
  c_ethanol_bg2_ta => carbon__mem % c_ethanol_bg2_ta % ptr ( ics : ice , iblk&
   )
  c_nonsoluble_bg2_ta => carbon__mem % c_nonsoluble_bg2_ta % ptr ( ics : ice ,&
   iblk )
  c_humus_2_ta => carbon__mem % c_humus_2_ta % ptr ( ics : ice , iblk )
!$acc data &
!$acc present(c_water_ag2_ta(:),c_water_bg2_ta(:),c_water_bg1_ta(:) &
!$acc ,current_c_state_sum(ics:ice),c_nonsoluble_bg2_ta(:),c_woods_ta(:) &
!$acc ,c_green_ta(:),c_crop_harvest_ta(:),c_nonsoluble_ag2_ta(:) &
!$acc ,c_ethanol_bg2_ta(:),c_acid_ag1_ta(:),c_humus_1_ta(:) &
!$acc ,c_nonsoluble_bg1_ta(:),c_humus_2_ta(:),c_nonsoluble_ag1_ta(:) &
!$acc ,c_ethanol_ag2_ta(:),c_acid_bg1_ta(:),c_ethanol_bg1_ta(:) &
!$acc ,c_reserve_ta(:),c_water_ag1_ta(:),c_acid_ag2_ta(:),c_ethanol_ag1_ta(:) &
!$acc ,c_acid_bg2_ta(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_2 = ics , ice
   current_c_state_sum ( claw_induction_2 ) = c_green_ta ( claw_induction_2 )&
    + c_woods_ta ( claw_induction_2 ) + c_reserve_ta ( claw_induction_2 ) +&
    c_crop_harvest_ta ( claw_induction_2 ) + c_acid_ag1_ta ( claw_induction_2&
    ) + c_water_ag1_ta ( claw_induction_2 ) + c_ethanol_ag1_ta (&
    claw_induction_2 ) + c_nonsoluble_ag1_ta ( claw_induction_2 ) +&
    c_acid_ag2_ta ( claw_induction_2 ) + c_water_ag2_ta ( claw_induction_2 ) +&
    c_ethanol_ag2_ta ( claw_induction_2 ) + c_nonsoluble_ag2_ta (&
    claw_induction_2 ) + c_acid_bg1_ta ( claw_induction_2 ) + c_water_bg1_ta (&
    claw_induction_2 ) + c_ethanol_bg1_ta ( claw_induction_2 ) +&
    c_nonsoluble_bg1_ta ( claw_induction_2 ) + c_acid_bg2_ta (&
    claw_induction_2 ) + c_water_bg2_ta ( claw_induction_2 ) +&
    c_ethanol_bg2_ta ( claw_induction_2 ) + c_nonsoluble_bg2_ta (&
    claw_induction_2 ) + c_humus_1_ta ( claw_induction_2 ) + c_humus_2_ta (&
    claw_induction_2 )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calculate_current_c_ta_state_sum

 SUBROUTINE calculate_current_c_ag_1_and_bg_sums ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , TARGET , INTENT(IN) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_carbon_interface" //&
  ":calculate_current_c_ag_1_and_bg_sums"
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: iblk
  REAL ( KIND= wp ) , POINTER :: c_bg_sum ( : )
  REAL ( KIND= wp ) , POINTER :: c_ag_sum_1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_2 ( : )
  INTEGER :: claw_induction_0

  SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
  TYPE IS ( t_carbon_memory )
   carbon__mem => mm
  END SELECT
  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  c_bg_sum => carbon__mem % c_bg_sum % ptr ( ics : ice , iblk )
  c_ag_sum_1 => carbon__mem % c_ag_sum_1 % ptr ( ics : ice , iblk )
  c_acid_ag1 => carbon__mem % c_acid_ag1 % ptr ( ics : ice , iblk )
  c_water_ag1 => carbon__mem % c_water_ag1 % ptr ( ics : ice , iblk )
  c_ethanol_ag1 => carbon__mem % c_ethanol_ag1 % ptr ( ics : ice , iblk )
  c_nonsoluble_ag1 => carbon__mem % c_nonsoluble_ag1 % ptr ( ics : ice , iblk&
   )
  c_acid_bg1 => carbon__mem % c_acid_bg1 % ptr ( ics : ice , iblk )
  c_water_bg1 => carbon__mem % c_water_bg1 % ptr ( ics : ice , iblk )
  c_ethanol_bg1 => carbon__mem % c_ethanol_bg1 % ptr ( ics : ice , iblk )
  c_nonsoluble_bg1 => carbon__mem % c_nonsoluble_bg1 % ptr ( ics : ice , iblk&
   )
  c_humus_1 => carbon__mem % c_humus_1 % ptr ( ics : ice , iblk )
  c_acid_bg2 => carbon__mem % c_acid_bg2 % ptr ( ics : ice , iblk )
  c_water_bg2 => carbon__mem % c_water_bg2 % ptr ( ics : ice , iblk )
  c_ethanol_bg2 => carbon__mem % c_ethanol_bg2 % ptr ( ics : ice , iblk )
  c_nonsoluble_bg2 => carbon__mem % c_nonsoluble_bg2 % ptr ( ics : ice , iblk&
   )
  c_humus_2 => carbon__mem % c_humus_2 % ptr ( ics : ice , iblk )
!$acc data &
!$acc present(c_bg_sum(ics:ice),c_acid_bg1(:),c_acid_bg2(:),c_ethanol_bg1(:) &
!$acc ,c_ethanol_bg2(:),c_ag_sum_1(ics:ice),c_nonsoluble_bg2(:) &
!$acc ,c_water_bg2(:),c_humus_2(:),c_water_bg1(:),c_humus_1(:),c_water_ag1(:) &
!$acc ,c_nonsoluble_bg1(:),c_acid_ag1(:),c_nonsoluble_ag1(:),c_ethanol_ag1(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_0 = ics , ice
   c_bg_sum ( claw_induction_0 ) = c_acid_bg1 ( claw_induction_0 ) +&
    c_water_bg1 ( claw_induction_0 ) + c_ethanol_bg1 ( claw_induction_0 ) +&
    c_nonsoluble_bg1 ( claw_induction_0 ) + c_acid_bg2 ( claw_induction_0 ) +&
    c_water_bg2 ( claw_induction_0 ) + c_ethanol_bg2 ( claw_induction_0 ) +&
    c_nonsoluble_bg2 ( claw_induction_0 ) + c_humus_1 ( claw_induction_0 ) +&
    c_humus_2 ( claw_induction_0 )
   c_ag_sum_1 ( claw_induction_0 ) = c_acid_ag1 ( claw_induction_0 ) +&
    c_water_ag1 ( claw_induction_0 ) + c_ethanol_ag1 ( claw_induction_0 ) +&
    c_nonsoluble_ag1 ( claw_induction_0 ) + c_humus_1 ( claw_induction_0 )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calculate_current_c_ag_1_and_bg_sums

 SUBROUTINE recalc_per_tile_vars ( tile , options , c_green , c_woods ,&
  c_reserve , c_crop_harvest , c_acid_ag1 , c_water_ag1 , c_ethanol_ag1 ,&
  c_nonsoluble_ag1 , c_acid_bg1 , c_water_bg1 , c_ethanol_bg1 ,&
  c_nonsoluble_bg1 , c_acid_ag2 , c_water_ag2 , c_ethanol_ag2 ,&
  c_nonsoluble_ag2 , c_acid_bg2 , c_water_bg2 , c_ethanol_bg2 ,&
  c_nonsoluble_bg2 , c_humus_1 , c_humus_2 , soil_respiration ,&
  npp_flux_correction , npp_pot_ydaymean , npp_act_ydaymean , gpp_ydaymean ,&
  root_exudates , cflux_c_green_2_herb , cflux_herb_2_atm ,&
  cflux_dist_greenreserve_2_soil , cflux_dist_woods_2_soil ,&
  cflux_fire_all_2_atm , cflux_c_greenwood_2_litter )
  USE mo_carbon_process , ONLY: get_per_tile

  CLASS ( t_jsb_tile_abstract ) , TARGET , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_green ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_woods ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_reserve ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_crop_harvest ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_acid_ag1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_water_ag1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_ethanol_ag1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_nonsoluble_ag1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_acid_bg1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_water_bg1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_ethanol_bg1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_nonsoluble_bg1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_acid_ag2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_water_ag2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_ethanol_ag2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_nonsoluble_ag2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_acid_bg2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_water_bg2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_ethanol_bg2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_nonsoluble_bg2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_humus_1 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: c_humus_2 ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: cflux_dist_greenreserve_2_soil&
   ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: cflux_dist_woods_2_soil ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: cflux_fire_all_2_atm ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: soil_respiration ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: npp_pot_ydaymean ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: npp_act_ydaymean ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: npp_flux_correction ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: gpp_ydaymean ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: cflux_c_greenwood_2_litter ( :&
   )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: root_exudates ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: cflux_c_green_2_herb ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: cflux_herb_2_atm ( : )
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_carbon_interface" //&
  ":recalc_per_tile_vars"
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: iblk
  REAL ( KIND= wp ) , POINTER :: c_green_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_reserve_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_woods_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_crop_harvest_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_1_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_2_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_dist_greenreserve_2_soil_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_dist_woods_2_soil_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_fire_all_2_atm_ta ( : )
  REAL ( KIND= wp ) , POINTER :: soil_respiration_ta ( : )
  REAL ( KIND= wp ) , POINTER :: npp_pot_ydaymean_ta ( : )
  REAL ( KIND= wp ) , POINTER :: npp_act_ydaymean_ta ( : )
  REAL ( KIND= wp ) , POINTER :: npp_flux_correction_ta ( : )
  REAL ( KIND= wp ) , POINTER :: gpp_ydaymean_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_c_greenwood_2_litter_ta ( : )
  REAL ( KIND= wp ) , POINTER :: root_exudates_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_c_green_2_herb_ta ( : )
  REAL ( KIND= wp ) , POINTER :: cflux_herb_2_atm_ta ( : )
  REAL ( KIND= wp ) , POINTER :: c_sum_natural_ta ( : )
  REAL ( KIND= wp ) , POINTER :: veg_fract_correction ( : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc_max ( : )

  SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
  TYPE IS ( t_carbon_memory )
   carbon__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  c_green_ta => carbon__mem % c_green_ta % ptr ( ics : ice , iblk )
  c_reserve_ta => carbon__mem % c_reserve_ta % ptr ( ics : ice , iblk )
  c_woods_ta => carbon__mem % c_woods_ta % ptr ( ics : ice , iblk )
  c_crop_harvest_ta => carbon__mem % c_crop_harvest_ta % ptr ( ics : ice ,&
   iblk )
  c_acid_ag1_ta => carbon__mem % c_acid_ag1_ta % ptr ( ics : ice , iblk )
  c_water_ag1_ta => carbon__mem % c_water_ag1_ta % ptr ( ics : ice , iblk )
  c_ethanol_ag1_ta => carbon__mem % c_ethanol_ag1_ta % ptr ( ics : ice , iblk&
   )
  c_nonsoluble_ag1_ta => carbon__mem % c_nonsoluble_ag1_ta % ptr ( ics : ice ,&
   iblk )
  c_acid_bg1_ta => carbon__mem % c_acid_bg1_ta % ptr ( ics : ice , iblk )
  c_water_bg1_ta => carbon__mem % c_water_bg1_ta % ptr ( ics : ice , iblk )
  c_ethanol_bg1_ta => carbon__mem % c_ethanol_bg1_ta % ptr ( ics : ice , iblk&
   )
  c_nonsoluble_bg1_ta => carbon__mem % c_nonsoluble_bg1_ta % ptr ( ics : ice ,&
   iblk )
  c_humus_1_ta => carbon__mem % c_humus_1_ta % ptr ( ics : ice , iblk )
  c_acid_ag2_ta => carbon__mem % c_acid_ag2_ta % ptr ( ics : ice , iblk )
  c_water_ag2_ta => carbon__mem % c_water_ag2_ta % ptr ( ics : ice , iblk )
  c_ethanol_ag2_ta => carbon__mem % c_ethanol_ag2_ta % ptr ( ics : ice , iblk&
   )
  c_nonsoluble_ag2_ta => carbon__mem % c_nonsoluble_ag2_ta % ptr ( ics : ice ,&
   iblk )
  c_acid_bg2_ta => carbon__mem % c_acid_bg2_ta % ptr ( ics : ice , iblk )
  c_water_bg2_ta => carbon__mem % c_water_bg2_ta % ptr ( ics : ice , iblk )
  c_ethanol_bg2_ta => carbon__mem % c_ethanol_bg2_ta % ptr ( ics : ice , iblk&
   )
  c_nonsoluble_bg2_ta => carbon__mem % c_nonsoluble_bg2_ta % ptr ( ics : ice ,&
   iblk )
  c_humus_2_ta => carbon__mem % c_humus_2_ta % ptr ( ics : ice , iblk )
  cflux_dist_greenreserve_2_soil_ta => carbon__mem %&
   cflux_dist_greenreserve_2_soil_ta % ptr ( ics : ice , iblk )
  cflux_dist_woods_2_soil_ta => carbon__mem % cflux_dist_woods_2_soil_ta % ptr&
   ( ics : ice , iblk )
  cflux_fire_all_2_atm_ta => carbon__mem % cflux_fire_all_2_atm_ta % ptr ( ics&
   : ice , iblk )
  soil_respiration_ta => carbon__mem % soil_respiration_ta % ptr ( ics : ice ,&
   iblk )
  npp_pot_ydaymean_ta => carbon__mem % npp_pot_ydaymean_ta % ptr ( ics : ice ,&
   iblk )
  npp_act_ydaymean_ta => carbon__mem % npp_act_ydaymean_ta % ptr ( ics : ice ,&
   iblk )
  npp_flux_correction_ta => carbon__mem % npp_flux_correction_ta % ptr ( ics :&
   ice , iblk )
  gpp_ydaymean_ta => carbon__mem % gpp_ydaymean_ta % ptr ( ics : ice , iblk )
  cflux_c_greenwood_2_litter_ta => carbon__mem % cflux_c_greenwood_2_litter_ta&
   % ptr ( ics : ice , iblk )
  root_exudates_ta => carbon__mem % root_exudates_ta % ptr ( ics : ice , iblk&
   )
  cflux_c_green_2_herb_ta => carbon__mem % cflux_c_green_2_herb_ta % ptr ( ics&
   : ice , iblk )
  cflux_herb_2_atm_ta => carbon__mem % cflux_herb_2_atm_ta % ptr ( ics : ice ,&
   iblk )
  c_sum_natural_ta => carbon__mem % c_sum_natural_ta % ptr ( ics : ice , iblk&
   )
  veg_fract_correction => pheno__mem % veg_fract_correction % ptr ( ics : ice&
   , iblk )
  fract_fpc_max => pheno__mem % fract_fpc_max % ptr ( ics : ice , iblk )
  IF ( present ( c_green ) ) THEN
   CALL get_per_tile ( c_green_ta ( : ) , c_green ( : ) , veg_fract_correction&
    ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_woods ) ) THEN
   CALL get_per_tile ( c_woods_ta ( : ) , c_woods ( : ) , veg_fract_correction&
    ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_reserve ) ) THEN
   CALL get_per_tile ( c_reserve_ta ( : ) , c_reserve ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_crop_harvest ) ) THEN
   CALL get_per_tile ( c_crop_harvest_ta ( : ) , c_crop_harvest ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_acid_ag1 ) ) THEN
   CALL get_per_tile ( c_acid_ag1_ta ( : ) , c_acid_ag1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_water_ag1 ) ) THEN
   CALL get_per_tile ( c_water_ag1_ta ( : ) , c_water_ag1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_ethanol_ag1 ) ) THEN
   CALL get_per_tile ( c_ethanol_ag1_ta ( : ) , c_ethanol_ag1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_nonsoluble_ag1 ) ) THEN
   CALL get_per_tile ( c_nonsoluble_ag1_ta ( : ) , c_nonsoluble_ag1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_acid_bg1 ) ) THEN
   CALL get_per_tile ( c_acid_bg1_ta ( : ) , c_acid_bg1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_water_bg1 ) ) THEN
   CALL get_per_tile ( c_water_bg1_ta ( : ) , c_water_bg1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_ethanol_bg1 ) ) THEN
   CALL get_per_tile ( c_ethanol_bg1_ta ( : ) , c_ethanol_bg1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_nonsoluble_bg1 ) ) THEN
   CALL get_per_tile ( c_nonsoluble_bg1_ta ( : ) , c_nonsoluble_bg1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_humus_1 ) ) THEN
   CALL get_per_tile ( c_humus_1_ta ( : ) , c_humus_1 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_acid_ag2 ) ) THEN
   CALL get_per_tile ( c_acid_ag2_ta ( : ) , c_acid_ag2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_water_ag2 ) ) THEN
   CALL get_per_tile ( c_water_ag2_ta ( : ) , c_water_ag2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_ethanol_ag2 ) ) THEN
   CALL get_per_tile ( c_ethanol_ag2_ta ( : ) , c_ethanol_ag2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_nonsoluble_ag2 ) ) THEN
   CALL get_per_tile ( c_nonsoluble_ag2_ta ( : ) , c_nonsoluble_ag2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_acid_bg2 ) ) THEN
   CALL get_per_tile ( c_acid_bg2_ta ( : ) , c_acid_bg2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_water_bg2 ) ) THEN
   CALL get_per_tile ( c_water_bg2_ta ( : ) , c_water_bg2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_ethanol_bg2 ) ) THEN
   CALL get_per_tile ( c_ethanol_bg2_ta ( : ) , c_ethanol_bg2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_nonsoluble_bg2 ) ) THEN
   CALL get_per_tile ( c_nonsoluble_bg2_ta ( : ) , c_nonsoluble_bg2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( c_humus_2 ) ) THEN
   CALL get_per_tile ( c_humus_2_ta ( : ) , c_humus_2 ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  CALL calculate_current_c_ta_state_sum ( tile , options , c_sum_natural_ta (&
   : ) )
  IF ( present ( cflux_dist_greenreserve_2_soil ) ) THEN
   CALL get_per_tile ( cflux_dist_greenreserve_2_soil_ta ( : ) ,&
    cflux_dist_greenreserve_2_soil ( : ) , veg_fract_correction ( : ) ,&
    fract_fpc_max ( : ) )
  END IF
  IF ( present ( cflux_dist_woods_2_soil ) ) THEN
   CALL get_per_tile ( cflux_dist_woods_2_soil_ta ( : ) ,&
    cflux_dist_woods_2_soil ( : ) , veg_fract_correction ( : ) , fract_fpc_max&
    ( : ) )
  END IF
  IF ( present ( cflux_fire_all_2_atm ) ) THEN
   CALL get_per_tile ( cflux_fire_all_2_atm_ta ( : ) , cflux_fire_all_2_atm (&
    : ) , veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( npp_pot_ydaymean ) ) THEN
   CALL get_per_tile ( npp_pot_ydaymean_ta ( : ) , npp_pot_ydaymean ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( npp_act_ydaymean ) ) THEN
   CALL get_per_tile ( npp_act_ydaymean_ta ( : ) , npp_act_ydaymean ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( soil_respiration ) ) THEN
   CALL get_per_tile ( soil_respiration_ta ( : ) , soil_respiration ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( npp_flux_correction ) ) THEN
   CALL get_per_tile ( npp_flux_correction_ta ( : ) , npp_flux_correction ( :&
    ) , veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( gpp_ydaymean ) ) THEN
   CALL get_per_tile ( gpp_ydaymean_ta ( : ) , gpp_ydaymean ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( cflux_c_greenwood_2_litter ) ) THEN
   CALL get_per_tile ( cflux_c_greenwood_2_litter_ta ( : ) ,&
    cflux_c_greenwood_2_litter ( : ) , veg_fract_correction ( : ) ,&
    fract_fpc_max ( : ) )
  END IF
  IF ( present ( root_exudates ) ) THEN
   CALL get_per_tile ( root_exudates_ta ( : ) , root_exudates ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( cflux_herb_2_atm ) ) THEN
   CALL get_per_tile ( cflux_herb_2_atm_ta ( : ) , cflux_herb_2_atm ( : ) ,&
    veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
  IF ( present ( cflux_c_green_2_herb ) ) THEN
   CALL get_per_tile ( cflux_c_green_2_herb_ta ( : ) , cflux_c_green_2_herb (&
    : ) , veg_fract_correction ( : ) , fract_fpc_max ( : ) )
  END IF
 END SUBROUTINE recalc_per_tile_vars

 SUBROUTINE check_carbon_conservation ( tile , options , old_c_state_sum ,&
  cflux , cconserve )

  CLASS ( t_jsb_tile_abstract ) , TARGET , INTENT(IN) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  REAL ( KIND= wp ) , INTENT(IN) :: old_c_state_sum ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: cflux ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: cconserve ( : )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_carbon_interface" //&
  ":check_carbon_conservation"
  INTEGER :: claw_induction_1

  CALL calculate_current_c_ta_state_sum ( tile , options , cconserve )
!$acc data present(cconserve(:),cflux(:),old_c_state_sum(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_1 = 1 , size ( cconserve , 1 )
   cconserve ( claw_induction_1 ) = cconserve ( claw_induction_1 ) -&
    old_c_state_sum ( claw_induction_1 ) - cflux ( claw_induction_1 )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE check_carbon_conservation

 SUBROUTINE aggregate_c_npp_pot_allocation ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_carbon_interface" //&
   ":aggregate_C_NPP_pot_allocation" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
  TYPE IS ( t_carbon_memory )
   carbon__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % lai_ydaymean , ics&
   , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_green_ta , ics ,&
   ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_woods_ta , ics ,&
   ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_reserve_ta , ics&
   , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_crop_harvest_ta&
   , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   c_litter_green_ag_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   c_litter_green_bg_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   c_litter_wood_ag_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   c_litter_wood_bg_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_acid_ag1_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_water_ag1_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_ethanol_ag1_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   c_nonsoluble_ag1_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_acid_bg1_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_water_bg1_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_ethanol_bg1_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   c_nonsoluble_bg1_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_humus_1_ta , ics&
   , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_acid_ag2_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_water_ag2_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_ethanol_ag2_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   c_nonsoluble_ag2_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_acid_bg2_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_water_bg2_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_ethanol_bg2_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   c_nonsoluble_bg2_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_humus_2_ta , ics&
   , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % c_sum_natural_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   cflux_c_greenwood_2_litter_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   soil_respiration_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   npp_flux_correction_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % excess_npp_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   npp_pot_ydaymean_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   npp_act_ydaymean_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % gpp_ydaymean_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % root_exudates_ta ,&
   ics , ice , iblk ,"mo_carbon_interface" //":aggregate_C_NPP_pot_allocation"&
   )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   cflux_c_green_2_herb_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   cconservation_calccpools , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem % cflux_npp_2_atm_ta&
   , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   cflux_soilresp_2_atm_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   cflux_herb_2_atm_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  CALL weighted_by_fract % aggregate ( tile , carbon__mem %&
   cflux_npp_2_atm_yday_ta , ics , ice , iblk ,"mo_carbon_interface" //&
  ":aggregate_C_NPP_pot_allocation" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_carbon_interface" //&
   ":aggregate_C_NPP_pot_allocation" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_c_npp_pot_allocation

 SUBROUTINE rescale_carbon_upon_reference_area_change ( tile , options ,&
  oldrefarea , newrefarea )
  USE mo_util , ONLY: real2string

  CLASS ( t_jsb_tile_abstract ) , TARGET , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  REAL ( KIND= wp ) , INTENT(IN) :: oldrefarea ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: newrefarea ( : )
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_carbon_interface" //&
  ":rescale_carbon_upon_reference_area_change"
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: iblk
  REAL ( KIND= wp ) , POINTER :: c_green ( : )
  REAL ( KIND= wp ) , POINTER :: c_reserve ( : )
  REAL ( KIND= wp ) , POINTER :: c_woods ( : )
  REAL ( KIND= wp ) , POINTER :: c_crop_harvest ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_ag2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_ag2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_ag2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_ag2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_acid_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_water_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_ethanol_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_nonsoluble_bg2 ( : )
  REAL ( KIND= wp ) , POINTER :: c_humus_2 ( : )

  IF ( any ( oldrefarea .le. 0 ) .OR. any ( oldrefarea .gt. 1 ) ) THEN
   CALL finish ( trim ("mo_carbon_interface" //&
   ":rescale_carbon_upon_reference_area_change" ) ,&
   "Violation of assertion: Reference areas (here old) need to be >0 and <=1 -&
   & min: " // real2string ( minval ( oldrefarea ) ) //"; max: " //&
    real2string ( maxval ( oldrefarea ) ) )
  END IF
  IF ( any ( newrefarea .le. 0 ) .OR. any ( newrefarea .gt. 1 ) ) THEN
   CALL finish ( trim ("mo_carbon_interface" //&
   ":rescale_carbon_upon_reference_area_change" ) ,&
   "Violation of assertion: Reference areas (here new) need to be >0 and <=1 -&
   & min: " // real2string ( minval ( newrefarea ) ) //"; max: " //&
    real2string ( maxval ( newrefarea ) ) )
  END IF
  SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
  TYPE IS ( t_carbon_memory )
   carbon__mem => mm
  END SELECT
  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  c_green => carbon__mem % c_green % ptr ( ics : ice , iblk )
  c_reserve => carbon__mem % c_reserve % ptr ( ics : ice , iblk )
  c_woods => carbon__mem % c_woods % ptr ( ics : ice , iblk )
  c_crop_harvest => carbon__mem % c_crop_harvest % ptr ( ics : ice , iblk )
  c_acid_ag1 => carbon__mem % c_acid_ag1 % ptr ( ics : ice , iblk )
  c_water_ag1 => carbon__mem % c_water_ag1 % ptr ( ics : ice , iblk )
  c_ethanol_ag1 => carbon__mem % c_ethanol_ag1 % ptr ( ics : ice , iblk )
  c_nonsoluble_ag1 => carbon__mem % c_nonsoluble_ag1 % ptr ( ics : ice , iblk&
   )
  c_acid_bg1 => carbon__mem % c_acid_bg1 % ptr ( ics : ice , iblk )
  c_water_bg1 => carbon__mem % c_water_bg1 % ptr ( ics : ice , iblk )
  c_ethanol_bg1 => carbon__mem % c_ethanol_bg1 % ptr ( ics : ice , iblk )
  c_nonsoluble_bg1 => carbon__mem % c_nonsoluble_bg1 % ptr ( ics : ice , iblk&
   )
  c_humus_1 => carbon__mem % c_humus_1 % ptr ( ics : ice , iblk )
  c_acid_ag2 => carbon__mem % c_acid_ag2 % ptr ( ics : ice , iblk )
  c_water_ag2 => carbon__mem % c_water_ag2 % ptr ( ics : ice , iblk )
  c_ethanol_ag2 => carbon__mem % c_ethanol_ag2 % ptr ( ics : ice , iblk )
  c_nonsoluble_ag2 => carbon__mem % c_nonsoluble_ag2 % ptr ( ics : ice , iblk&
   )
  c_acid_bg2 => carbon__mem % c_acid_bg2 % ptr ( ics : ice , iblk )
  c_water_bg2 => carbon__mem % c_water_bg2 % ptr ( ics : ice , iblk )
  c_ethanol_bg2 => carbon__mem % c_ethanol_bg2 % ptr ( ics : ice , iblk )
  c_nonsoluble_bg2 => carbon__mem % c_nonsoluble_bg2 % ptr ( ics : ice , iblk&
   )
  c_humus_2 => carbon__mem % c_humus_2 % ptr ( ics : ice , iblk )
  c_green = c_green * ( oldrefarea / newrefarea )
  c_woods = c_woods * ( oldrefarea / newrefarea )
  c_reserve = c_reserve * ( oldrefarea / newrefarea )
  c_crop_harvest = c_crop_harvest * ( oldrefarea / newrefarea )
  c_acid_ag1 = c_acid_ag1 * ( oldrefarea / newrefarea )
  c_water_ag1 = c_water_ag1 * ( oldrefarea / newrefarea )
  c_ethanol_ag1 = c_ethanol_ag1 * ( oldrefarea / newrefarea )
  c_nonsoluble_ag1 = c_nonsoluble_ag1 * ( oldrefarea / newrefarea )
  c_acid_bg1 = c_acid_bg1 * ( oldrefarea / newrefarea )
  c_water_bg1 = c_water_bg1 * ( oldrefarea / newrefarea )
  c_ethanol_bg1 = c_ethanol_bg1 * ( oldrefarea / newrefarea )
  c_nonsoluble_bg1 = c_nonsoluble_bg1 * ( oldrefarea / newrefarea )
  c_humus_1 = c_humus_1 * ( oldrefarea / newrefarea )
  c_acid_ag2 = c_acid_ag2 * ( oldrefarea / newrefarea )
  c_water_ag2 = c_water_ag2 * ( oldrefarea / newrefarea )
  c_ethanol_ag2 = c_ethanol_ag2 * ( oldrefarea / newrefarea )
  c_nonsoluble_ag2 = c_nonsoluble_ag2 * ( oldrefarea / newrefarea )
  c_acid_bg2 = c_acid_bg2 * ( oldrefarea / newrefarea )
  c_water_bg2 = c_water_bg2 * ( oldrefarea / newrefarea )
  c_ethanol_bg2 = c_ethanol_bg2 * ( oldrefarea / newrefarea )
  c_nonsoluble_bg2 = c_nonsoluble_bg2 * ( oldrefarea / newrefarea )
  c_humus_2 = c_humus_2 * ( oldrefarea / newrefarea )
  CALL calculate_current_c_ag_1_and_bg_sums ( tile , options )
  CALL recalc_per_tile_vars ( tile , options , c_green = c_green ( : ) ,&
   c_woods = c_woods ( : ) , c_reserve = c_reserve ( : ) , c_crop_harvest =&
   c_crop_harvest ( : ) , c_acid_ag1 = c_acid_ag1 ( : ) , c_water_ag1 =&
   c_water_ag1 ( : ) , c_ethanol_ag1 = c_ethanol_ag1 ( : ) , c_nonsoluble_ag1&
   = c_nonsoluble_ag1 ( : ) , c_acid_ag2 = c_acid_ag2 ( : ) , c_water_ag2 =&
   c_water_ag2 ( : ) , c_ethanol_ag2 = c_ethanol_ag2 ( : ) , c_nonsoluble_ag2&
   = c_nonsoluble_ag2 ( : ) , c_acid_bg1 = c_acid_bg1 ( : ) , c_water_bg1 =&
   c_water_bg1 ( : ) , c_ethanol_bg1 = c_ethanol_bg1 ( : ) , c_nonsoluble_bg1&
   = c_nonsoluble_bg1 ( : ) , c_acid_bg2 = c_acid_bg2 ( : ) , c_water_bg2 =&
   c_water_bg2 ( : ) , c_ethanol_bg2 = c_ethanol_bg2 ( : ) , c_nonsoluble_bg2&
   = c_nonsoluble_bg2 ( : ) , c_humus_1 = c_humus_1 ( : ) , c_humus_2 =&
   c_humus_2 ( : ) )
 END SUBROUTINE rescale_carbon_upon_reference_area_change

 SUBROUTINE yday_carbon_conservation_test ( tile )
  USE mo_carbon_constants , ONLY: molarmassco2_kg , sec_per_day
  USE mo_jsb_grid , ONLY: get_grid

  CLASS ( t_jsb_tile_abstract ) , TARGET , INTENT(IN) :: tile
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  TYPE ( t_l2a_memory ) , POINTER :: l2a__mem
  REAL ( KIND= wp ) , POINTER :: carbon_conservation_test ( : )
  REAL ( KIND= wp ) , POINTER :: yday_c_state_sum ( : )
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_jsb_grid ) , POINTER :: grid
  TYPE ( t_jsb_task_options ) :: options
  REAL ( KIND= wp ) , ALLOCATABLE :: c_flux_yday ( : )
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: iblk
  INTEGER :: nc

  SELECT TYPE ( mm => tile % mem ( l2a_ ) % p )
  TYPE IS ( t_l2a_memory )
   l2a__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
  TYPE IS ( t_carbon_memory )
   carbon__mem => mm
  END SELECT
  model => get_model ( tile % owner_model_id )
  grid => get_grid ( model % grid_id )
  options % ics = 1
  options % ice = grid % nproma
  ics = options % ics
  ice = options % ice
  nc = ice - ics + 1
  ALLOCATE ( c_flux_yday ( nc ) )
  DO iblk = 1 , grid % nblks , 1
   CALL model % set_options ( iblk = iblk )
   options % iblk = iblk
   IF ( iblk == grid % nblks ) THEN
    options % ice = grid % npromz
    ice = grid % npromz
   END IF
   carbon_conservation_test => l2a__mem % carbon_conservation_test % ptr ( ics&
    : ice , iblk )
   yday_c_state_sum => l2a__mem % yday_c_state_sum % ptr ( ics : ice , iblk )
   IF ( ( .NOT. all ( carbon_conservation_test .eq. ( - 999.0_wp ) ) ) ) THEN
    c_flux_yday ( ics : ice ) = ( - ( carbon__mem % cflux_npp_2_atm_yday_ta %&
     ptr ( ics : ice , iblk ) + carbon__mem % cflux_soilresp_2_atm_ta % ptr (&
     ics : ice , iblk ) + carbon__mem % cflux_herb_2_atm_ta % ptr ( ics : ice&
     , iblk ) ) / molarmassco2_kg * sec_per_day )
    IF ( model % processes ( disturb_ ) % p % config % active ) THEN
     c_flux_yday ( ics : ice ) = c_flux_yday ( ics : ice ) - carbon__mem %&
      cflux_fire_all_2_atm_ta % ptr ( ics : ice , iblk )
    END IF
    CALL calculate_current_c_ta_state_sum ( tile , options ,&
     carbon_conservation_test )
    carbon_conservation_test = carbon_conservation_test - yday_c_state_sum -&
     c_flux_yday ( ics : ice )
   ELSE
    carbon_conservation_test = 0._wp
   END IF
   CALL calculate_current_c_ta_state_sum ( tile , options , yday_c_state_sum )
  END DO
  DEALLOCATE ( c_flux_yday )
 END SUBROUTINE yday_carbon_conservation_test

END MODULE mo_carbon_interface

