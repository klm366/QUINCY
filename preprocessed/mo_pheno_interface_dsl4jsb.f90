MODULE mo_pheno_interface
 USE mo_jsb_control , ONLY: debug_on
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message , finish
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract , t_jsb_aggregator
 USE mo_jsb_process_class , ONLY: t_jsb_process
 USE mo_jsb_task_class , ONLY: t_jsb_process_task , t_jsb_task_options
 USE mo_jsb_process_class , ONLY: seb_ , a2l_ , hydro_ , assimi_ , pheno_ ,&
  carbon_
 USE mo_pheno_config_class , ONLY: t_pheno_config
 USE mo_a2l_memory_class , ONLY: t_a2l_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 USE mo_pheno_memory_class , ONLY: t_pheno_memory
 USE mo_assimi_memory_class , ONLY: t_assimi_memory
 USE mo_carbon_memory_class , ONLY: t_carbon_memory
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_pheno_interface"
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_phenology
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_phenology
  PROCEDURE , NOPASS :: aggregate => aggregate_phenology
 END TYPE tsk_phenology
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_foliage_projected_cover
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_fpc
  PROCEDURE , NOPASS :: aggregate => aggregate_fpc
 END TYPE tsk_foliage_projected_cover
 PRIVATE :: create_task_phenology
 PRIVATE :: create_task_foliage_projected_cover
 PRIVATE :: update_phenology
 PRIVATE :: update_phenology_logrop
 PRIVATE :: update_phenology_climatology
 PRIVATE :: aggregate_phenology
 PRIVATE :: update_fpc
 PRIVATE :: aggregate_fpc
 PUBLIC :: register_pheno_tasks
 INTERFACE tsk_phenology
  PROCEDURE create_task_phenology
 END INTERFACE
 INTERFACE tsk_foliage_projected_cover
  PROCEDURE create_task_foliage_projected_cover
 END INTERFACE

CONTAINS
 FUNCTION create_task_phenology ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_phenology :: return_ptr )
  CALL return_ptr % construct ( name ="phenology" , process_id = pheno_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_phenology

 FUNCTION create_task_foliage_projected_cover ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_foliage_projected_cover :: return_ptr )
  CALL return_ptr % construct ( name ="foliage_projected_cover" , process_id =&
   pheno_ , owner_model_id = model_id )
 END FUNCTION create_task_foliage_projected_cover

 SUBROUTINE register_pheno_tasks ( this , model_id )

  CLASS ( t_jsb_process ) , INTENT(INOUT) :: this
  INTEGER , INTENT(IN) :: model_id

  CALL this % register_task ( tsk_phenology ( model_id ) )
  CALL this % register_task ( tsk_foliage_projected_cover ( model_id ) )
 END SUBROUTINE register_pheno_tasks

 SUBROUTINE update_phenology ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_pheno_config ) , POINTER :: pheno__conf
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_interface" //&
  ":update_phenology"

  IF ( ( .NOT. tile % is_process_active ( pheno_ ) ) ) THEN
   RETURN
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( pheno_ ) % p % config )
  TYPE IS ( t_pheno_config )
   pheno__conf => cc
  END SELECT
  IF ( tile % lcts ( 1 ) % lib_id == 0 .AND. trim ( pheno__conf % scheme ) /=&
  "climatology" ) THEN
   CALL message ( trim ("mo_pheno_interface" //":update_phenology" ) ,&
   "Vegetation tile has no PFTs, using climatology for LAI" )
   CALL update_phenology_climatology ( tile , options )
  ELSE
   SELECT CASE ( trim ( pheno__conf % scheme ) )
   CASE ("logrop" )
    CALL update_phenology_logrop ( tile , options )
   CASE ("climatology" )
    CALL update_phenology_climatology ( tile , options )
   CASE DEFAULT
    CALL finish ("Invalid scheme for phenology" )
   END SELECT
  END IF
 END SUBROUTINE update_phenology

 SUBROUTINE update_phenology_climatology ( tile , options )
  USE mo_jsb_time , ONLY: get_time_interpolation_weights

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: wgt1
  REAL ( KIND= wp ) :: wgt2
  INTEGER :: nmw1
  INTEGER :: nmw2
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: lai_mon ( : , : )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_interface" //&
  ":update_phenology_climatology"

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( pheno_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":update_phenology_climatology"&
    ) ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  lai => pheno__mem % lai % ptr ( ics : ice , iblk )
  lai_mon ( 1 : , 0 : ) => pheno__mem % lai_mon ( ics : ice , iblk , : )
  CALL get_time_interpolation_weights ( wgt1 , wgt2 , nmw1 , nmw2 )
  lai ( : ) = real ( wgt1 , wp ) * lai_mon ( : , nmw1 ) + real ( wgt2 , wp ) *&
   lai_mon ( : , nmw2 )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":update_phenology_climatology"&
    ) ,"Finished." )
  END IF
 END SUBROUTINE update_phenology_climatology

 SUBROUTINE update_phenology_logrop ( tile , options )
  USE mo_carbon_interface , ONLY: rescale_carbon_upon_reference_area_change ,&
   calculate_current_c_ta_state_sum , check_carbon_conservation
  USE mo_pheno_process , ONLY: calc_summergreen_phenology ,&
   calc_evergreen_phenology , calc_crop_phenology , calc_raingreen_phenology ,&
   calc_grass_phenology , track_max_green_pool ,&
   derive_maxlai_according_to_allometric_relationships
  USE mo_jsb_time_iface , ONLY: get_time_current , get_year_day
  USE mo_jsb_time , ONLY: is_newyear , is_newday , timestep_in_days
  USE mo_pheno_constants , ONLY: phenoparam
  USE mo_jsb_grid_class , ONLY: t_jsb_grid
  USE mo_jsb_grid , ONLY: get_grid
  USE mo_physical_constants , ONLY: tmelt

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: day_of_year
  REAL ( KIND= wp ) :: w_soil_filling ( 1 : options % nc )
  REAL ( KIND= wp ) :: t_air_in_celcius ( 1 : options % nc )
  REAL ( KIND= wp ) :: laimax_corrected ( 1 : options % nc )
  REAL ( KIND= wp ) :: veg_fract_correction_old ( 1 : options % nc )
  REAL ( KIND= wp ) :: old_c_state_sum ( 1 : options % nc )
  REAL ( KIND= wp ) :: dummy_flux ( 1 : options % nc )
  TYPE ( t_jsb_grid ) , POINTER :: grid
  TYPE ( t_pheno_config ) , POINTER :: pheno__conf
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  REAL ( KIND= wp ) , POINTER :: t_air ( : )
  REAL ( KIND= wp ) , POINTER :: npp_pot_rate_ca ( : )
  REAL ( KIND= wp ) , POINTER :: w_soil_rel ( : )
  REAL ( KIND= wp ) , POINTER :: lat ( : )
  REAL ( KIND= wp ) , POINTER :: heat_sum_sg ( : )
  REAL ( KIND= wp ) , POINTER :: heat_sum_eg ( : )
  REAL ( KIND= wp ) , POINTER :: heat_sum_crp ( : )
  REAL ( KIND= wp ) , POINTER :: heat_sum_winter ( : )
  REAL ( KIND= wp ) , POINTER :: days_since_growth_begin_sg ( : )
  REAL ( KIND= wp ) , POINTER :: days_since_growth_begin_eg ( : )
  REAL ( KIND= wp ) , POINTER :: growth_phase_sg ( : )
  REAL ( KIND= wp ) , POINTER :: growth_phase_eg ( : )
  REAL ( KIND= wp ) , POINTER :: growth_phase_crp ( : )
  REAL ( KIND= wp ) , POINTER :: previous_day_temp_mean ( : )
  REAL ( KIND= wp ) , POINTER :: pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , POINTER :: chill_days_sg ( : )
  REAL ( KIND= wp ) , POINTER :: chill_days_eg ( : )
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: day_npp_sum ( : )
  REAL ( KIND= wp ) , POINTER :: previous_day_temp_min ( : )
  REAL ( KIND= wp ) , POINTER :: previous_day_temp_max ( : )
  REAL ( KIND= wp ) , POINTER :: previous_day_npp_pot_rate_ca ( : )
  REAL ( KIND= wp ) , POINTER :: c_green ( : )
  REAL ( KIND= wp ) , POINTER :: c_reserve ( : )
  REAL ( KIND= wp ) , POINTER :: c_woods ( : )
  REAL ( KIND= wp ) , POINTER :: current_max_green ( : )
  REAL ( KIND= wp ) , POINTER :: veg_carbon_at_max_green ( : )
  REAL ( KIND= wp ) , POINTER :: number_of_individuals ( : )
  REAL ( KIND= wp ) , POINTER :: biomass_per_individual ( : )
  REAL ( KIND= wp ) , POINTER :: maxlai_allom ( : )
  REAL ( KIND= wp ) , POINTER :: veg_fract_correction ( : )
  REAL ( KIND= wp ) , POINTER :: cconservation_allom ( : )
  LOGICAL :: l_newday
  REAL ( KIND= wp ) :: timestep
  INTEGER :: i
  REAL ( KIND= wp ) :: maxlai
  REAL ( KIND= wp ) :: lai_negligible_tmp
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_interface" //&
  ":update_phenology_logrop"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_5
  INTEGER :: claw_induction_6
  INTEGER :: claw_induction_8
  INTEGER :: claw_induction_9

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( pheno_ ) ) ) THEN
   RETURN
  END IF
!$acc data create(w_soil_filling, t_air_in_celcius) &
!$acc create(laimax_corrected, veg_fract_correction_old) &
!$acc create(old_c_state_sum, dummy_flux)
  IF ( ( .NOT. tile % is_process_active ( assimi_ ) ) ) THEN
   CALL finish ("mo_pheno_interface" //":update_phenology_logrop" ,&
   "The LOGROP phenology model depends on the assimilation process being activ&
   &ated, too" )
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":update_phenology_logrop" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( pheno_ ) % p % config )
  TYPE IS ( t_pheno_config )
   pheno__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  grid => get_grid ( model % grid_id )
  lat => grid % lat ( ics : ice , iblk )
  IF ( pheno__conf % l_forestregrowth ) THEN
   IF ( ( .NOT. tile % is_process_active ( carbon_ ) ) ) THEN
    CALL finish ( trim ("mo_pheno_interface" //":update_phenology_logrop" ) ,&
    "Forest regrowth can only be switched on when the CARBON process is active&
    &, too" )
   END IF
   SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
   TYPE IS ( t_carbon_memory )
    carbon__mem => mm
   END SELECT
   c_green => carbon__mem % c_green % ptr ( ics : ice , iblk )
   c_reserve => carbon__mem % c_reserve % ptr ( ics : ice , iblk )
   c_woods => carbon__mem % c_woods % ptr ( ics : ice , iblk )
   current_max_green => carbon__mem % current_max_green % ptr ( ics : ice ,&
    iblk )
   veg_carbon_at_max_green => carbon__mem % veg_carbon_at_max_green % ptr (&
    ics : ice , iblk )
   number_of_individuals => pheno__mem % number_of_individuals % ptr ( ics :&
    ice , iblk )
   biomass_per_individual => pheno__mem % biomass_per_individual % ptr ( ics :&
    ice , iblk )
   maxlai_allom => pheno__mem % maxlai_allom % ptr ( ics : ice , iblk )
   veg_fract_correction => pheno__mem % veg_fract_correction % ptr ( ics : ice&
    , iblk )
   cconservation_allom => pheno__mem % cconservation_allom % ptr ( ics : ice ,&
    iblk )
  END IF
  t_air => a2l__mem % t_air % ptr ( ics : ice , iblk )
  npp_pot_rate_ca => assimi__mem % npp_pot_rate_ca % ptr ( ics : ice , iblk )
  w_soil_rel => hydro__mem % w_soil_rel % ptr ( ics : ice , iblk )
  previous_day_temp_mean => seb__mem % previous_day_temp_mean % ptr ( ics :&
   ice , iblk )
  previous_day_temp_min => seb__mem % previous_day_temp_min % ptr ( ics : ice&
   , iblk )
  previous_day_temp_max => seb__mem % previous_day_temp_max % ptr ( ics : ice&
   , iblk )
  pseudo_soil_temp => seb__mem % pseudo_soil_temp % ptr ( ics : ice , iblk )
  lai => pheno__mem % lai % ptr ( ics : ice , iblk )
  heat_sum_sg => pheno__mem % heat_sum_sg % ptr ( ics : ice , iblk )
  heat_sum_eg => pheno__mem % heat_sum_eg % ptr ( ics : ice , iblk )
  heat_sum_crp => pheno__mem % heat_sum_crp % ptr ( ics : ice , iblk )
  heat_sum_winter => pheno__mem % heat_sum_winter % ptr ( ics : ice , iblk )
  days_since_growth_begin_sg => pheno__mem % days_since_growth_begin_sg % ptr&
   ( ics : ice , iblk )
  days_since_growth_begin_eg => pheno__mem % days_since_growth_begin_eg % ptr&
   ( ics : ice , iblk )
  growth_phase_sg => pheno__mem % growth_phase_sg % ptr ( ics : ice , iblk )
  growth_phase_eg => pheno__mem % growth_phase_eg % ptr ( ics : ice , iblk )
  growth_phase_crp => pheno__mem % growth_phase_crp % ptr ( ics : ice , iblk )
  chill_days_sg => pheno__mem % chill_days_sg % ptr ( ics : ice , iblk )
  chill_days_eg => pheno__mem % chill_days_eg % ptr ( ics : ice , iblk )
  day_npp_sum => pheno__mem % day_npp_sum % ptr ( ics : ice , iblk )
  previous_day_npp_pot_rate_ca => pheno__mem % previous_day_npp_pot_rate_ca %&
   ptr ( ics : ice , iblk )
  IF ( ( .NOT. tile % is_vegetation ) ) THEN
   CALL finish ( trim ("mo_pheno_interface" //":update_phenology_logrop" ) ,&
   "The programm tries to calculate phenology on a non vegetated tile..." )
  END IF
  l_newday = is_newday ( )
  timestep = timestep_in_days ( )
!$acc data present(t_air(:),t_air_in_celcius(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_2 = 1 , size ( t_air_in_celcius , 1 )
   t_air_in_celcius ( claw_induction_2 ) = t_air ( claw_induction_2 ) - tmelt
  END DO
!$acc end parallel
!$acc end data
  day_of_year = int ( get_year_day ( get_time_current ( ) ) )
!$acc data present(w_soil_rel(:),w_soil_filling(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_1 = 1 , size ( w_soil_filling , 1 )
   w_soil_filling ( claw_induction_1 ) = w_soil_rel ( claw_induction_1 )
  END DO
!$acc end parallel
!$acc end data
  IF ( l_newday ) THEN
!$acc data present(day_npp_sum(:),previous_day_npp_pot_rate_ca(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_0 = 1 , size ( previous_day_npp_pot_rate_ca , 1 )
    previous_day_npp_pot_rate_ca ( claw_induction_0 ) = day_npp_sum (&
     claw_induction_0 ) / 86400._wp
    day_npp_sum ( claw_induction_0 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( pheno__conf % l_forestregrowth .AND. model % lctlib ( tile % lcts ( 1 )&
   % lib_id ) % forestflag ) THEN
   IF ( is_newyear ( ) ) THEN
!$acc data present(veg_fract_correction_old(:),veg_fract_correction(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_6 = 1 , size ( veg_fract_correction_old , 1 )
     veg_fract_correction_old ( claw_induction_6 ) = veg_fract_correction (&
      claw_induction_6 )
    END DO
!$acc end parallel
!$acc end data
    CALL derive_maxlai_according_to_allometric_relationships ( model % lctlib&
     ( tile % lcts ( 1 ) % lib_id ) % alpha_nr_ind , model % lctlib ( tile %&
     lcts ( 1 ) % lib_id ) % beta_nr_ind , model % lctlib ( tile % lcts ( 1 )&
     % lib_id ) % alpha_leaf , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
     beta_leaf , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
     clumpinessfactor , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
     specificleafarea_c , c_woods ( : ) , current_max_green ( : ) ,&
     veg_carbon_at_max_green ( : ) , number_of_individuals ( : ) ,&
     biomass_per_individual ( : ) , veg_fract_correction ( : ) , maxlai_allom&
     ( : ) )
    CALL calculate_current_c_ta_state_sum ( tile , options , old_c_state_sum (&
     : ) )
    CALL rescale_carbon_upon_reference_area_change ( tile , options ,&
     veg_fract_correction_old , veg_fract_correction )
!$acc data present(dummy_flux(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_5 = 1 , size ( dummy_flux , 1 )
     dummy_flux ( claw_induction_5 ) = 0._wp
    END DO
!$acc end parallel
!$acc end data
    CALL check_carbon_conservation ( tile , options , old_c_state_sum ( : ) ,&
     dummy_flux ( : ) , cconservation_allom ( : ) )
   END IF
   IF ( is_newday ( ) ) THEN
    CALL track_max_green_pool ( c_green ( : ) , c_reserve ( : ) , c_woods ( :&
     ) , current_max_green ( : ) , veg_carbon_at_max_green ( : ) )
   END IF
  END IF
  IF ( pheno__conf % l_forestregrowth ) THEN
!$acc data present(pheno__mem%maxlai_allom%ptr(ics:ice,iblk),laimax_corrected(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_9 = 1 , size ( laimax_corrected , 1 )
    laimax_corrected ( claw_induction_9 ) = pheno__mem % maxlai_allom % ptr (&
     claw_induction_9 , iblk )
   END DO
!$acc end parallel
!$acc end data
  ELSE
   maxlai = model % lctlib ( tile % lcts ( 1 ) % lib_id ) % maxlai
!$acc data present(laimax_corrected(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_8 = 1 , size ( laimax_corrected , 1 )
    laimax_corrected ( claw_induction_8 ) = maxlai
   END DO
!$acc end parallel
!$acc end data
  END IF
  lai_negligible_tmp = phenoparam % all % lai_negligible
!$acc data present(laimax_corrected, lai)
!$acc parallel loop gang vector
  DO i = 1 , size ( laimax_corrected ) , 1
   IF ( laimax_corrected ( i ) < lai_negligible_tmp ) THEN
    laimax_corrected ( i ) = lai_negligible_tmp
    lai ( i ) = lai_negligible_tmp - epsilon ( 1.0_wp )
   END IF
  END DO
!$acc end data
  IF ( model % lctlib ( tile % lcts ( 1 ) % lib_id ) % phenologytype == 2 )&
   THEN
   IF ( debug_on ( ) .AND. iblk == 1 ) THEN
    CALL message ("     " ,"... summergreen ..." )
   END IF
   CALL calc_summergreen_phenology ( day_of_year , l_newday , timestep , lat (&
    : ) , previous_day_temp_mean ( : ) , pseudo_soil_temp ( : ) ,&
    laimax_corrected ( : ) , days_since_growth_begin_sg ( : ) , chill_days_sg&
    ( : ) , heat_sum_sg ( : ) , growth_phase_sg ( : ) , lai ( : ) )
  END IF
  IF ( model % lctlib ( tile % lcts ( 1 ) % lib_id ) % phenologytype == 1 )&
   THEN
   IF ( debug_on ( ) .AND. iblk == 1 ) THEN
    CALL message ("     " ,"... evergreen ..." )
   END IF
   CALL calc_evergreen_phenology ( day_of_year , l_newday , timestep , lat ( :&
    ) , pseudo_soil_temp ( : ) , laimax_corrected ( : ) ,&
    days_since_growth_begin_eg ( : ) , chill_days_eg ( : ) , heat_sum_eg ( : )&
    , growth_phase_eg ( : ) , lai ( : ) )
  END IF
  IF ( model % lctlib ( tile % lcts ( 1 ) % lib_id ) % phenologytype == 5 )&
   THEN
   IF ( debug_on ( ) .AND. iblk == 1 ) THEN
    CALL message ("     " ,"... crop ..." )
   END IF
   CALL calc_crop_phenology ( day_of_year , l_newday , timestep ,&
    npp_pot_rate_ca ( : ) , lat ( : ) , pseudo_soil_temp ( : ) ,&
    laimax_corrected ( : ) , previous_day_npp_pot_rate_ca ( : ) ,&
    w_soil_filling ( : ) , previous_day_temp_min ( : ) , previous_day_temp_max&
    ( : ) , model % lctlib ( tile % lcts ( 1 ) % lib_id ) % specificleafarea_c&
    , heat_sum_crp ( : ) , heat_sum_winter , growth_phase_crp ( : ) , lai ( :&
    ) )
  END IF
  IF ( model % lctlib ( tile % lcts ( 1 ) % lib_id ) % phenologytype == 3 )&
   THEN
   IF ( debug_on ( ) .AND. iblk == 1 ) THEN
    CALL message ("     " ,"... raingreen ..." )
   END IF
   CALL calc_raingreen_phenology ( timestep , w_soil_filling ( : ) ,&
    laimax_corrected ( : ) , previous_day_npp_pot_rate_ca ( : ) , lai ( : ) )
  END IF
  IF ( model % lctlib ( tile % lcts ( 1 ) % lib_id ) % phenologytype == 4 )&
   THEN
   IF ( debug_on ( ) .AND. iblk == 1 ) THEN
    CALL message ("     " ,"... grasses ..." )
   END IF
   CALL calc_grass_phenology ( t_air_in_celcius ( : ) , timestep ,&
    w_soil_filling ( : ) , laimax_corrected ( : ) ,&
    previous_day_npp_pot_rate_ca ( : ) , lai ( : ) )
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":update_phenology_logrop" ) ,&
   "Finished." )
  END IF
!$acc end data
 END SUBROUTINE update_phenology_logrop

 SUBROUTINE aggregate_phenology ( tile , options )
  USE mo_jsb_time , ONLY: is_newyear

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_pheno_config ) , POINTER :: pheno__conf
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_interface" //&
  ":aggregate_phenology"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":aggregate_phenology" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( pheno_ ) % p % config )
  TYPE IS ( t_pheno_config )
   pheno__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , pheno__mem % lai , ics , ice ,&
   iblk ,"mo_pheno_interface" //":aggregate_phenology" )
  IF ( pheno__conf % l_forestregrowth ) THEN
   IF ( is_newyear ( ) ) THEN
    CALL weighted_by_fract % aggregate ( tile , pheno__mem % maxlai_allom ,&
     ics , ice , iblk ,"mo_pheno_interface" //":aggregate_phenology" )
    CALL weighted_by_fract % aggregate ( tile , pheno__mem %&
     cconservation_allom , ics , ice , iblk ,"mo_pheno_interface" //&
    ":aggregate_phenology" )
   END IF
  END IF
  IF ( trim ( pheno__conf % scheme ) =="logrop" ) THEN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":aggregate_phenology" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_phenology

 SUBROUTINE update_fpc ( tile , options )
  USE mo_pheno_process , ONLY: get_foliage_projected_cover
  USE mo_jsb_time , ONLY: get_time_interpolation_weights

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_pheno_config ) , POINTER :: pheno__conf
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc ( : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc_max ( : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc_mon ( : , : )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: wgt1
  REAL ( KIND= wp ) :: wgt2
  REAL ( KIND= wp ) :: clumpinessfactor
  INTEGER :: nmw1
  INTEGER :: nmw2
  TYPE ( t_jsb_model ) , POINTER :: model
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_interface" //&
  ":update_fpc"
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4
  INTEGER :: claw_induction_7

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( pheno_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":update_fpc" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( pheno_ ) % p % config )
  TYPE IS ( t_pheno_config )
   pheno__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  lai => pheno__mem % lai % ptr ( ics : ice , iblk )
  fract_fpc_max => pheno__mem % fract_fpc_max % ptr ( ics : ice , iblk )
  fract_fpc => pheno__mem % fract_fpc % ptr ( ics : ice , iblk )
  IF ( trim ( pheno__conf % scheme ) =="climatology" .OR. tile % lcts ( 1 ) %&
   lib_id == 0 ) THEN
   IF ( model % config % l_compat401 ) THEN
!$acc data present(lai(:),fract_fpc_max(:),fract_fpc(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_4 = 1 , size ( fract_fpc , 1 )
     fract_fpc ( claw_induction_4 ) = get_foliage_projected_cover (&
      fract_fpc_max ( claw_induction_4 ) , lai ( claw_induction_4 ) , 1._wp )
    END DO
!$acc end parallel
!$acc end data
   ELSE
    fract_fpc_mon ( 1 : , 0 : ) => pheno__mem % fract_fpc_mon ( ics : ice ,&
     iblk , : )
    CALL get_time_interpolation_weights ( wgt1 , wgt2 , nmw1 , nmw2 )
!$acc data present(fract_fpc_mon(:,nmw2),fract_fpc_mon(:,nmw1),fract_fpc(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_3 = 1 , size ( fract_fpc , 1 )
     fract_fpc ( claw_induction_3 ) = real ( wgt1 , wp ) * fract_fpc_mon (&
      claw_induction_3 , nmw1 ) + real ( wgt2 , wp ) * fract_fpc_mon (&
      claw_induction_3 , nmw2 )
    END DO
!$acc end parallel
!$acc end data
   END IF
  ELSE
   clumpinessfactor = model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    clumpinessfactor
!$acc data present(lai(:),fract_fpc_max(:),fract_fpc(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_7 = 1 , size ( fract_fpc , 1 )
    fract_fpc ( claw_induction_7 ) = get_foliage_projected_cover (&
     fract_fpc_max ( claw_induction_7 ) , lai ( claw_induction_7 ) ,&
     clumpinessfactor )
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":update_fpc" ) ,"Finished." )
  END IF
 END SUBROUTINE update_fpc

 SUBROUTINE aggregate_fpc ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_interface" //&
  ":aggregate_fpc"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":aggregate_fpc" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , pheno__mem % fract_fpc , ics ,&
   ice , iblk ,"mo_pheno_interface" //":aggregate_fpc" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_pheno_interface" //":aggregate_fpc" ) ,"Finished."&
    )
  END IF
 END SUBROUTINE aggregate_fpc

END MODULE mo_pheno_interface

