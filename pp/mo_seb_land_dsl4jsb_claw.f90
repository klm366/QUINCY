MODULE mo_seb_land
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract
 USE mo_jsb_task_class , ONLY: t_jsb_task_options
 USE mo_jsb_lct_class , ONLY: glacier_type , contains_lct
 USE mo_jsb_control , ONLY: debug_on
 USE mo_jsb_time , ONLY: is_time_experiment_start
 USE mo_jsb_process_class , ONLY: seb_ , rad_ , turb_ , a2l_ , sse_ , hydro_
 USE mo_sse_config_class , ONLY: t_sse_config
 USE mo_a2l_memory_class , ONLY: t_a2l_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_rad_memory_class , ONLY: t_rad_memory
 USE mo_turb_memory_class , ONLY: t_turb_memory
 USE mo_sse_memory_class , ONLY: t_sse_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 PRIVATE :: surface_temp_implicit
 PRIVATE :: calc_previous_day_variables
 PRIVATE :: calc_pseudo_soil_temp
 PUBLIC :: update_surface_energy_land
 PUBLIC :: update_asselin_land
 PUBLIC :: update_surface_fluxes_land
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_seb_land"

CONTAINS
 SUBROUTINE update_surface_energy_land ( tile , options )
  USE mo_sse_process , ONLY: get_liquid_max
  USE mo_phy_schemes , ONLY: qsat , q_effective , surface_dry_static_energy ,&
   heat_transfer_coef
  USE mo_jsb_physical_constants , ONLY: tmelt , rhoh2o , tpfac2 , tpfac3
  USE mo_sse_constants , ONLY: snow_depth_min
  USE mo_jsb_time , ONLY: is_newday , timesteps_per_day

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  TYPE ( t_turb_memory ) , POINTER :: turb__mem
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  LOGICAL :: lstart
  REAL ( KIND= wp ) , POINTER :: t_air ( : )
  REAL ( KIND= wp ) , POINTER :: day_temp_min ( : )
  REAL ( KIND= wp ) , POINTER :: day_temp_max ( : )
  REAL ( KIND= wp ) , POINTER :: day_temp_sum ( : )
  REAL ( KIND= wp ) , POINTER :: previous_day_temp_mean ( : )
  REAL ( KIND= wp ) , POINTER :: previous_day_temp_min ( : )
  REAL ( KIND= wp ) , POINTER :: previous_day_temp_max ( : )
  REAL ( KIND= wp ) , POINTER :: f_pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , POINTER :: n_pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , POINTER :: pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , POINTER :: t ( : )
  REAL ( KIND= wp ) , POINTER :: t_old ( : )
  REAL ( KIND= wp ) , POINTER :: t_unfilt ( : )
  REAL ( KIND= wp ) , POINTER :: t_unfilt_old ( : )
  REAL ( KIND= wp ) , POINTER :: t_eff ( : )
  REAL ( KIND= wp ) , POINTER :: qsat_star ( : )
  REAL ( KIND= wp ) , POINTER :: s_star ( : )
  REAL ( KIND= wp ) , POINTER :: heat_cap ( : )
  REAL ( KIND= wp ) , POINTER :: forc_hflx ( : )
  REAL ( KIND= wp ) , POINTER :: press_srf ( : )
  REAL ( KIND= wp ) , POINTER :: q_air ( : )
  REAL ( KIND= wp ) , POINTER :: fact_qsat_srf ( : )
  REAL ( KIND= wp ) , POINTER :: fact_q_air ( : )
  REAL ( KIND= wp ) , POINTER :: drag_srf ( : )
  REAL ( KIND= wp ) , POINTER :: rad_srf_net ( : )
  REAL ( KIND= wp ) , POINTER :: grnd_hflx ( : )
  REAL ( KIND= wp ) , POINTER :: hcap_grnd ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow ( : )
  REAL ( KIND= wp ) , POINTER :: snow_soil_dens ( : )
  REAL ( KIND= wp ) , POINTER :: soil_depth_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_ice_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_porosity_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: matrix_pot_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: bclapp_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: t_acoef ( : )
  REAL ( KIND= wp ) , POINTER :: t_bcoef ( : )
  REAL ( KIND= wp ) , POINTER :: q_acoef ( : )
  REAL ( KIND= wp ) , POINTER :: q_bcoef ( : )
  REAL ( KIND= wp ) :: fract ( 1 : options % nc )
  REAL ( KIND= wp ) :: qsat_srf_old ( 1 : options % nc )
  REAL ( KIND= wp ) :: s_old ( 1 : options % nc )
  REAL ( KIND= wp ) :: dqdt ( 1 : options % nc )
  REAL ( KIND= wp ) :: t2s_conv ( 1 : options % nc )
  REAL ( KIND= wp ) :: heat_tcoef ( 1 : options % nc )
  REAL ( KIND= wp ) :: t_star ( 1 : options % nc )
  REAL ( KIND= wp ) :: liquid_max ( 1 : options % nc )
  REAL ( KIND= wp ) :: t_air_in_celcius ( 1 : options % nc )
  LOGICAL :: is_glacier ( 1 : options % nc )
  LOGICAL :: has_snow_layer ( 1 : options % nc )
  REAL ( KIND= wp ) :: soil_depth_sl_vol_porosity_sl_tmp ( 1 : size (&
   liquid_max , dim = 1 ) )
  REAL ( KIND= wp ) :: w_soil_critical_tmp
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: nc
  INTEGER :: i
  REAL ( KIND= wp ) :: steplen
  INTEGER :: ilct
  TYPE ( t_jsb_model ) , POINTER :: model
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_land" //&
  ":update_surface_energy_land"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4
  INTEGER :: claw_induction_9
  INTEGER :: claw_induction_10
  INTEGER :: claw_induction_11
  INTEGER :: claw_induction_12
  INTEGER :: claw_induction_13
  INTEGER :: claw_induction_15
  INTEGER :: claw_induction_16
  INTEGER :: claw_induction_17
  INTEGER :: claw_induction_18
  INTEGER :: claw_induction_19

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  nc = options % nc
  steplen = options % steplen
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_land" //":update_surface_energy_land" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( sse_ ) % p % config )
  TYPE IS ( t_sse_config )
   sse__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( turb_ ) % p )
  TYPE IS ( t_turb_memory )
   turb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  t_air => a2l__mem % t_air % ptr ( ics : ice , iblk )
  t_acoef => a2l__mem % t_acoef % ptr ( ics : ice , iblk )
  t_bcoef => a2l__mem % t_bcoef % ptr ( ics : ice , iblk )
  q_acoef => a2l__mem % q_acoef % ptr ( ics : ice , iblk )
  q_bcoef => a2l__mem % q_bcoef % ptr ( ics : ice , iblk )
  press_srf => a2l__mem % press_srf % ptr ( ics : ice , iblk )
  q_air => a2l__mem % q_air % ptr ( ics : ice , iblk )
  drag_srf => a2l__mem % drag_srf % ptr ( ics : ice , iblk )
  previous_day_temp_mean => seb__mem % previous_day_temp_mean % ptr ( ics :&
   ice , iblk )
  day_temp_sum => seb__mem % day_temp_sum % ptr ( ics : ice , iblk )
  day_temp_min => seb__mem % day_temp_min % ptr ( ics : ice , iblk )
  day_temp_max => seb__mem % day_temp_max % ptr ( ics : ice , iblk )
  previous_day_temp_min => seb__mem % previous_day_temp_min % ptr ( ics : ice&
   , iblk )
  previous_day_temp_max => seb__mem % previous_day_temp_max % ptr ( ics : ice&
   , iblk )
  n_pseudo_soil_temp => seb__mem % n_pseudo_soil_temp % ptr ( ics : ice , iblk&
   )
  f_pseudo_soil_temp => seb__mem % f_pseudo_soil_temp % ptr ( ics : ice , iblk&
   )
  pseudo_soil_temp => seb__mem % pseudo_soil_temp % ptr ( ics : ice , iblk )
  t => seb__mem % t % ptr ( ics : ice , iblk )
  t_old => seb__mem % t_old % ptr ( ics : ice , iblk )
  t_unfilt => seb__mem % t_unfilt % ptr ( ics : ice , iblk )
  t_unfilt_old => seb__mem % t_unfilt_old % ptr ( ics : ice , iblk )
  t_eff => seb__mem % t_eff % ptr ( ics : ice , iblk )
  qsat_star => seb__mem % qsat_star % ptr ( ics : ice , iblk )
  s_star => seb__mem % s_star % ptr ( ics : ice , iblk )
  heat_cap => seb__mem % heat_cap % ptr ( ics : ice , iblk )
  forc_hflx => seb__mem % forc_hflx % ptr ( ics : ice , iblk )
  fact_q_air => turb__mem % fact_q_air % ptr ( ics : ice , iblk )
  fact_qsat_srf => turb__mem % fact_qsat_srf % ptr ( ics : ice , iblk )
  rad_srf_net => rad__mem % rad_srf_net % ptr ( ics : ice , iblk )
  w_snow => hydro__mem % w_snow % ptr ( ics : ice , iblk )
  fract_snow => hydro__mem % fract_snow % ptr ( ics : ice , iblk )
  snow_soil_dens => hydro__mem % snow_soil_dens % ptr ( ics : ice , iblk )
  soil_depth_sl => hydro__mem % soil_depth_sl % ptr ( ics : ice , : , iblk )
  w_soil_sl => hydro__mem % w_soil_sl % ptr ( ics : ice , : , iblk )
  w_ice_sl => hydro__mem % w_ice_sl % ptr ( ics : ice , : , iblk )
  vol_porosity_sl => hydro__mem % vol_porosity_sl % ptr ( ics : ice , : , iblk&
   )
  matrix_pot_sl => hydro__mem % matrix_pot_sl % ptr ( ics : ice , : , iblk )
  bclapp_sl => hydro__mem % bclapp_sl % ptr ( ics : ice , : , iblk )
  grnd_hflx => sse__mem % grnd_hflx % ptr ( ics : ice , iblk )
  hcap_grnd => sse__mem % hcap_grnd % ptr ( ics : ice , iblk )
!$acc data create(fract,qsat_srf_old,s_old,dqdt,t2s_conv,heat_tcoef, &
!$acc t_star,liquid_max,t_air_in_celcius,is_glacier,has_snow_layer, &
!$acc soil_depth_sl_vol_porosity_sl_tmp)
!$acc data present(is_glacier)
  IF ( contains_lct ( tile % lcts , glacier_type ) ) THEN
   DO ilct = 1 , size ( tile % lcts ) , 1
    IF ( tile % lcts ( ilct ) % id == glacier_type ) THEN
     is_glacier ( : ) = tile % lcts ( ilct ) % fract ( ics : ice , iblk ) >&
      0._wp
     EXIT
    END IF
   END DO
  ELSE
   is_glacier ( : ) = .FALSE.
  END IF
!$acc update device(is_glacier) !
!$acc end data
!$acc data &
!$acc present(soil_depth_sl_vol_porosity_sl_tmp(:),vol_porosity_sl(:,1) &
!$acc ,soil_depth_sl(:,1))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_19 = 1 , size ( soil_depth_sl_vol_porosity_sl_tmp , 1 )
   soil_depth_sl_vol_porosity_sl_tmp ( claw_induction_19 ) = soil_depth_sl (&
    claw_induction_19 , 1 ) * vol_porosity_sl ( claw_induction_19 , 1 )
  END DO
!$acc end parallel
!$acc end data
!$acc data present(t_air(:),t_air_in_celcius(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_18 = 1 , size ( t_air_in_celcius , 1 )
   t_air_in_celcius ( claw_induction_18 ) = t_air ( claw_induction_18 ) -&
    tmelt
  END DO
!$acc end parallel
!$acc end data
  CALL calc_previous_day_variables ( is_newday ( options % datetime , options&
   % dtime ) , t_air_in_celcius ( : ) , timesteps_per_day ( options % dtime )&
   , previous_day_temp_mean ( : ) , day_temp_sum ( : ) , previous_day_temp_min&
   ( : ) , day_temp_min ( : ) , previous_day_temp_max ( : ) , day_temp_max ( :&
   ) )
  CALL calc_pseudo_soil_temp ( t_air_in_celcius ( : ) , n_pseudo_soil_temp ( :&
   ) , f_pseudo_soil_temp ( : ) , pseudo_soil_temp ( : ) )
!$acc data present(t_unfilt_old(:),t_unfilt(:),t(:),t_old(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_0 = 1 , size ( t_unfilt_old , 1 )
   t_unfilt_old ( claw_induction_0 ) = t_unfilt ( claw_induction_0 )
   t_old ( claw_induction_0 ) = t ( claw_induction_0 )
  END DO
!$acc end parallel
!$acc end data
!$acc data &
!$acc present(forc_hflx(:),grnd_hflx(:),qsat_srf_old(:),heat_cap(:) &
!$acc ,hcap_grnd(:),t_old(:),press_srf(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_1 = 1 , size ( qsat_srf_old , 1 )
   qsat_srf_old ( claw_induction_1 ) = qsat ( t_old ( claw_induction_1 ) ,&
    press_srf ( claw_induction_1 ) )
   heat_cap ( claw_induction_1 ) = hcap_grnd ( claw_induction_1 )
   forc_hflx ( claw_induction_1 ) = grnd_hflx ( claw_induction_1 )
  END DO
!$acc end parallel
!$acc end data
  IF ( is_time_experiment_start ( options % datetime ) ) THEN
   lstart = .TRUE.
  ELSE
   lstart = .FALSE.
  END IF
!$acc data &
!$acc present(dqdt(:),qsat_srf_old(:),t(:),fact_q_air(:),t_old(:) &
!$acc ,heat_tcoef(:),fact_qsat_srf(:),press_srf(:),q_air(:),s_old(:) &
!$acc ,drag_srf(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_2 = 1 , size ( dqdt , 1 )
   dqdt ( claw_induction_2 ) = ( qsat ( t ( claw_induction_2 ) + 0.001_wp ,&
    press_srf ( claw_induction_2 ) ) - qsat_srf_old ( claw_induction_2 ) ) *&
    1000._wp
   s_old ( claw_induction_2 ) = surface_dry_static_energy ( t_old (&
    claw_induction_2 ) , q_effective ( qsat_srf_old ( claw_induction_2 ) ,&
    q_air ( claw_induction_2 ) , fact_qsat_srf ( claw_induction_2 ) ,&
    fact_q_air ( claw_induction_2 ) ) )
   heat_tcoef ( claw_induction_2 ) = heat_transfer_coef ( drag_srf (&
    claw_induction_2 ) , steplen )
  END DO
!$acc end parallel
!$acc end data
!$acc data present(t2s_conv(:),t_old(:),s_old(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_17 = 1 , size ( t2s_conv , 1 )
   t2s_conv ( claw_induction_17 ) = s_old ( claw_induction_17 ) / t_old (&
    claw_induction_17 )
  END DO
!$acc end parallel
!$acc end data
  IF ( lstart ) THEN
!$acc data present(s_old(:),s_star(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_13 = 1 , size ( s_star , 1 )
    s_star ( claw_induction_13 ) = s_old ( claw_induction_13 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
   CALL surface_temp_implicit ( nc , steplen , t2s_conv ( : ) , t_acoef ( : )&
    , t_bcoef ( : ) , q_acoef ( : ) , q_bcoef ( : ) , s_old ( : ) ,&
    qsat_srf_old ( : ) , dqdt ( : ) , rad_srf_net ( : ) , forc_hflx ( : ) ,&
    heat_tcoef ( : ) , fact_q_air ( : ) , fact_qsat_srf ( : ) , fract_snow ( :&
    ) , heat_cap ( : ) , s_star ( : ) )
  END IF
!$acc data present(fract)
  CALL tile % parent_tile % get_fraction ( ics , ice , iblk , fract = fract (&
   : ) )
!$acc update device(fract)
!$acc end data
!$acc data present(fract,s_star,t2s_conv)
!$acc parallel loop gang vector
  DO i = 1 , size ( s_star , dim = 1 ) , 1
   IF ( fract ( i ) < epsilon ( 1._wp ) ) THEN
    s_star ( i ) = 280._wp * t2s_conv ( i )
   END IF
  END DO
!$acc end data
!$acc data present(t_unfilt(:),t2s_conv(:),t_star(:),s_old(:),s_star(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_3 = 1 , size ( t_unfilt , 1 )
   t_unfilt ( claw_induction_3 ) = ( tpfac2 * s_star ( claw_induction_3 ) +&
    tpfac3 * s_old ( claw_induction_3 ) ) / t2s_conv ( claw_induction_3 )
   t_star ( claw_induction_3 ) = s_star ( claw_induction_3 ) / t2s_conv (&
    claw_induction_3 )
  END DO
!$acc end parallel
!$acc end data
  w_soil_critical_tmp = sse__conf % w_soil_critical
  IF ( sse__conf % l_snow ) THEN
!$acc data present(snow_soil_dens(:),w_snow(:),has_snow_layer(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_12 = 1 , size ( has_snow_layer , 1 )
    has_snow_layer ( claw_induction_12 ) = w_snow ( claw_induction_12 ) *&
     rhoh2o / snow_soil_dens ( claw_induction_12 ) > snow_depth_min
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(w_snow(:),has_snow_layer(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_11 = 1 , size ( has_snow_layer , 1 )
    has_snow_layer ( claw_induction_11 ) = w_snow ( claw_induction_11 ) >&
     w_soil_critical_tmp
   END DO
!$acc end parallel
!$acc end data
  END IF
!$acc data present(has_snow_layer,is_glacier,t_star)
!$acc parallel loop gang vector
  DO i = 1 , size ( t_star , dim = 1 ) , 1
   IF ( has_snow_layer ( i ) .OR. is_glacier ( i ) ) THEN
    t_star ( i ) = min ( t_star ( i ) , tmelt )
   END IF
  END DO
!$acc end data
  IF ( sse__conf % l_freeze ) THEN
   IF ( sse__conf % l_supercool ) THEN
!$acc data &
!$acc present(soil_depth_sl_vol_porosity_sl_tmp(:),bclapp_sl(:,1),t_star(:) &
!$acc ,matrix_pot_sl(:,1),liquid_max(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_4 = 1 , size ( liquid_max , 1 )
     liquid_max ( claw_induction_4 ) = get_liquid_max ( t_star (&
      claw_induction_4 ) , soil_depth_sl_vol_porosity_sl_tmp (&
      claw_induction_4 ) , matrix_pot_sl ( claw_induction_4 , 1 ) , bclapp_sl&
      ( claw_induction_4 , 1 ) )
    END DO
!$acc end parallel
!$acc end data
   ELSE
    liquid_max ( : ) = 0._wp
   END IF
!$acc data present(has_snow_layer,is_glacier,w_soil_sl,liquid_max,t_star,w_ice_sl)
!$acc parallel loop gang vector
   DO i = 1 , size ( t_star , dim = 1 ) , 1
    IF ( ( .NOT. has_snow_layer ( i ) ) .AND. ( .NOT. is_glacier ( i ) ) )&
     THEN
     IF ( w_soil_sl ( i , 1 ) > max ( liquid_max ( i ) , w_soil_critical_tmp )&
      ) THEN
      t_star ( i ) = max ( t_star ( i ) , tmelt )
     END IF
     IF ( w_ice_sl ( i , 1 ) > w_soil_critical_tmp ) THEN
      t_star ( i ) = min ( t_star ( i ) , tmelt )
     END IF
    END IF
   END DO
!$acc end data
  END IF
!$acc data present(t2s_conv(:),t_star(:),s_star(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_16 = 1 , size ( s_star , 1 )
   s_star ( claw_induction_16 ) = t_star ( claw_induction_16 ) * t2s_conv (&
    claw_induction_16 )
  END DO
!$acc end parallel
!$acc end data
  IF ( lstart ) THEN
!$acc data present(qsat_srf_old(:),qsat_star(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_10 = 1 , size ( qsat_star , 1 )
    qsat_star ( claw_induction_10 ) = qsat_srf_old ( claw_induction_10 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data &
!$acc present(dqdt(:),t2s_conv(:),qsat_srf_old(:),s_old(:),qsat_star(:) &
!$acc ,s_star(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_9 = 1 , size ( qsat_star , 1 )
    qsat_star ( claw_induction_9 ) = qsat_srf_old ( claw_induction_9 ) + dqdt&
     ( claw_induction_9 ) * ( s_star ( claw_induction_9 ) - s_old (&
     claw_induction_9 ) ) / t2s_conv ( claw_induction_9 )
   END DO
!$acc end parallel
!$acc end data
  END IF
!$acc data present(t_eff(:),t_star(:),t_old(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_15 = 1 , size ( t_eff , 1 )
   t_eff ( claw_induction_15 ) = ( t_old ( claw_induction_15 ) ** ( 3 ) * (&
    4._wp * t_star ( claw_induction_15 ) - 3._wp * t_old ( claw_induction_15 )&
    ) ) ** ( 0.25_wp )
  END DO
!$acc end parallel
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_land" //":update_surface_energy_land" ) ,&
   "Finished." )
  END IF
!$acc end data
 END SUBROUTINE update_surface_energy_land

 SUBROUTINE calc_previous_day_variables ( is_newday , t_air_in_celcius ,&
  time_steps_per_day , previous_day_temp_mean , day_temp_sum ,&
  previous_day_temp_min , day_temp_min , previous_day_temp_max , day_temp_max&
  )

  LOGICAL , INTENT(IN) :: is_newday
  REAL ( KIND= wp ) , INTENT(IN) :: t_air_in_celcius ( : )
  INTEGER , INTENT(IN) :: time_steps_per_day
  REAL ( KIND= wp ) , INTENT(INOUT) :: previous_day_temp_mean ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: previous_day_temp_min ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: previous_day_temp_max ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: day_temp_sum ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: day_temp_min ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: day_temp_max ( : )
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_air_in_celcius,previous_day_temp_mean,previous_day_temp_min &
!$acc ,previous_day_temp_max,day_temp_sum,day_temp_min,day_temp_max)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( previous_day_temp_max , 1 ) , 1
   IF ( is_newday ) THEN
    previous_day_temp_mean ( jsb_hori ) = day_temp_sum ( jsb_hori ) /&
     time_steps_per_day
    day_temp_sum ( jsb_hori ) = t_air_in_celcius ( jsb_hori )
    previous_day_temp_min ( jsb_hori ) = day_temp_min ( jsb_hori )
    day_temp_min ( jsb_hori ) = t_air_in_celcius ( jsb_hori )
    previous_day_temp_max ( jsb_hori ) = day_temp_max ( jsb_hori )
    day_temp_max ( jsb_hori ) = t_air_in_celcius ( jsb_hori )
   ELSE
    day_temp_sum ( jsb_hori ) = day_temp_sum ( jsb_hori ) + t_air_in_celcius (&
     jsb_hori )
    IF ( day_temp_min ( jsb_hori ) > t_air_in_celcius ( jsb_hori ) ) THEN
     day_temp_min ( jsb_hori ) = t_air_in_celcius ( jsb_hori )
    END IF
    IF ( day_temp_max ( jsb_hori ) < t_air_in_celcius ( jsb_hori ) ) THEN
     day_temp_max ( jsb_hori ) = t_air_in_celcius ( jsb_hori )
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_previous_day_variables

 SUBROUTINE calc_pseudo_soil_temp ( t_air_in_celcius , n_pseudo_soil_temp ,&
  f_pseudo_soil_temp , pseudo_soil_temp )

  REAL ( KIND= wp ) , INTENT(IN) :: t_air_in_celcius ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: n_pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: f_pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: pseudo_soil_temp ( : )
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_air_in_celcius,n_pseudo_soil_temp,f_pseudo_soil_temp &
!$acc ,pseudo_soil_temp)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( n_pseudo_soil_temp , 1 ) , 1
   pseudo_soil_temp ( jsb_hori ) = t_air_in_celcius ( jsb_hori ) /&
    n_pseudo_soil_temp ( jsb_hori ) + f_pseudo_soil_temp ( jsb_hori ) *&
    pseudo_soil_temp ( jsb_hori )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_pseudo_soil_temp

 SUBROUTINE surface_temp_implicit ( klon , steplen , pcp , pescoe , pfscoe ,&
  peqcoe , pfqcoe , psold , pqsold , pdqsold , pnetrad , pgrdfl , pcfh , pcair&
  , pcsat , pfracsu , pgrdcap , psnew )
  USE mo_jsb_physical_constants , ONLY: zemiss_def , stbo , alv , als ,&
   cvdifts

  INTEGER , INTENT(IN) :: klon
  REAL ( KIND= wp ) , INTENT(IN) :: steplen
  REAL ( KIND= wp ) , INTENT(IN) :: pcp ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pfscoe ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pescoe ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pfqcoe ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: peqcoe ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: psold ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pqsold ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pdqsold ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pnetrad ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pgrdfl ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pcfh ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pcair ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pcsat ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pfracsu ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(IN) :: pgrdcap ( 1 : klon )
  REAL ( KIND= wp ) , INTENT(OUT) :: psnew ( 1 : klon )
  REAL ( KIND= wp ) :: zcolin ( 1 : klon )
  REAL ( KIND= wp ) :: zcohfl ( 1 : klon )
  REAL ( KIND= wp ) :: zcoind ( 1 : klon )
  REAL ( KIND= wp ) :: zicp ( 1 : klon )
  REAL ( KIND= wp ) :: zca ( 1 : klon )
  REAL ( KIND= wp ) :: zcs ( 1 : klon )
  REAL ( KIND= wp ) :: pdt
  INTEGER :: i
  INTEGER :: claw_induction_5

!$acc data &
!$acc create(zcolin,zcohfl,zcoind,zicp,zca,zcs)
  pdt = cvdifts * steplen
!$acc data &
!$acc present(zcohfl(:),zcolin(:),psold(:),pqsold(:),pgrdfl(:),pcair(:) &
!$acc ,pnetrad(:),zicp(:),zca(:),psnew(:),pcsat(:),pescoe(:),zcoind(:) &
!$acc ,pcfh(:),pfscoe(:),pfracsu(:),pgrdcap(:),zcs(:),pfqcoe(:),peqcoe(:) &
!$acc ,pdqsold(:),pcp(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_5 = 1 , size ( zicp , 1 )
   zicp ( claw_induction_5 ) = 1._wp / pcp ( claw_induction_5 )
   zca ( claw_induction_5 ) = als * pfracsu ( claw_induction_5 ) + alv * (&
    pcair ( claw_induction_5 ) - pfracsu ( claw_induction_5 ) )
   zcs ( claw_induction_5 ) = als * pfracsu ( claw_induction_5 ) + alv * (&
    pcsat ( claw_induction_5 ) - pfracsu ( claw_induction_5 ) )
   zcolin ( claw_induction_5 ) = pgrdcap ( claw_induction_5 ) * zicp (&
    claw_induction_5 ) + pdt * ( zicp ( claw_induction_5 ) * 4._wp *&
    zemiss_def * stbo * ( zicp ( claw_induction_5 ) * psold ( claw_induction_5&
    ) ) ** ( 3 ) - pcfh ( claw_induction_5 ) * ( zca ( claw_induction_5 ) *&
    peqcoe ( claw_induction_5 ) - zcs ( claw_induction_5 ) ) * zicp (&
    claw_induction_5 ) * pdqsold ( claw_induction_5 ) )
   zcohfl ( claw_induction_5 ) = ( - pdt * pcfh ( claw_induction_5 ) * (&
    pescoe ( claw_induction_5 ) - 1._wp ) )
   zcoind ( claw_induction_5 ) = pdt * ( pnetrad ( claw_induction_5 ) + pcfh (&
    claw_induction_5 ) * pfscoe ( claw_induction_5 ) + pcfh ( claw_induction_5&
    ) * ( ( zca ( claw_induction_5 ) * peqcoe ( claw_induction_5 ) - zcs (&
    claw_induction_5 ) ) * pqsold ( claw_induction_5 ) + zca (&
    claw_induction_5 ) * pfqcoe ( claw_induction_5 ) ) + pgrdfl (&
    claw_induction_5 ) )
   psnew ( claw_induction_5 ) = ( zcolin ( claw_induction_5 ) * psold (&
    claw_induction_5 ) + zcoind ( claw_induction_5 ) ) / ( zcolin (&
    claw_induction_5 ) + zcohfl ( claw_induction_5 ) )
  END DO
!$acc end parallel
!$acc end data
!$acc end data
 END SUBROUTINE surface_temp_implicit

 SUBROUTINE update_asselin_land ( tile , options )
  USE mo_jsb_time , ONLY: get_asselin_coef

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  REAL ( KIND= wp ) , POINTER :: t ( : )
  REAL ( KIND= wp ) , POINTER :: t_filt ( : )
  REAL ( KIND= wp ) , POINTER :: t_old ( : )
  REAL ( KIND= wp ) , POINTER :: t_unfilt ( : )
  REAL ( KIND= wp ) , POINTER :: t_unfilt_old ( : )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: eps
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_land" //&
  ":update_asselin_land"
  INTEGER :: claw_induction_7
  INTEGER :: claw_induction_8
  INTEGER :: claw_induction_14

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_land" //":update_asselin_land" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  t => seb__mem % t % ptr ( ics : ice , iblk )
  t_filt => seb__mem % t_filt % ptr ( ics : ice , iblk )
  t_old => seb__mem % t_old % ptr ( ics : ice , iblk )
  t_unfilt => seb__mem % t_unfilt % ptr ( ics : ice , iblk )
  t_unfilt_old => seb__mem % t_unfilt_old % ptr ( ics : ice , iblk )
  eps = get_asselin_coef ( )
  IF ( eps > 0._wp ) THEN
!$acc data present(t_unfilt_old(:),t_filt(:),t_unfilt(:),t_old(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_8 = 1 , size ( t_filt , 1 )
    t_filt ( claw_induction_8 ) = t_unfilt_old ( claw_induction_8 ) + eps * (&
     t_old ( claw_induction_8 ) - 2._wp * t_unfilt_old ( claw_induction_8 ) +&
     t_unfilt ( claw_induction_8 ) )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(t_filt(:),t_unfilt(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_7 = 1 , size ( t_filt , 1 )
    t_filt ( claw_induction_7 ) = t_unfilt ( claw_induction_7 )
   END DO
!$acc end parallel
!$acc end data
  END IF
!$acc data present(t_filt(:),t(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_14 = 1 , size ( t , 1 )
   t ( claw_induction_14 ) = t_filt ( claw_induction_14 )
  END DO
!$acc end parallel
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_land" //":update_asselin_land" ) ,"Finished."&
    )
  END IF
 END SUBROUTINE update_asselin_land

 SUBROUTINE update_surface_fluxes_land ( tile , options )
  USE mo_phy_schemes , ONLY: heat_transfer_coef
  USE mo_jsb_physical_constants , ONLY: alv , als

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  REAL ( KIND= wp ) , POINTER :: s_star ( : )
  REAL ( KIND= wp ) , POINTER :: t_acoef ( : )
  REAL ( KIND= wp ) , POINTER :: t_bcoef ( : )
  REAL ( KIND= wp ) , POINTER :: drag_srf ( : )
  REAL ( KIND= wp ) , POINTER :: latent_hflx ( : )
  REAL ( KIND= wp ) , POINTER :: sensible_hflx ( : )
  REAL ( KIND= wp ) , POINTER :: latent_hflx_lnd ( : )
  REAL ( KIND= wp ) , POINTER :: sensible_hflx_lnd ( : )
  REAL ( KIND= wp ) , POINTER :: evapopot ( : )
  REAL ( KIND= wp ) , POINTER :: evapotrans ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow ( : )
  REAL ( KIND= wp ) , ALLOCATABLE :: s_air ( : )
  REAL ( KIND= wp ) , ALLOCATABLE :: heat_tcoef ( : )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: nc
  REAL ( KIND= wp ) :: steplen
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_land" //&
  ":update_surface_fluxes_land"
  INTEGER :: claw_induction_6

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  nc = options % nc
  steplen = options % steplen
  IF ( ( .NOT. tile % is_process_active ( seb_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_land" //":update_surface_fluxes_land" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  t_acoef => a2l__mem % t_acoef % ptr ( ics : ice , iblk )
  t_bcoef => a2l__mem % t_bcoef % ptr ( ics : ice , iblk )
  drag_srf => a2l__mem % drag_srf % ptr ( ics : ice , iblk )
  s_star => seb__mem % s_star % ptr ( ics : ice , iblk )
  evapotrans => hydro__mem % evapotrans % ptr ( ics : ice , iblk )
  evapopot => hydro__mem % evapopot % ptr ( ics : ice , iblk )
  fract_snow => hydro__mem % fract_snow % ptr ( ics : ice , iblk )
  sensible_hflx => seb__mem % sensible_hflx % ptr ( ics : ice , iblk )
  latent_hflx => seb__mem % latent_hflx % ptr ( ics : ice , iblk )
  sensible_hflx_lnd => seb__mem % sensible_hflx_lnd % ptr ( ics : ice , iblk )
  latent_hflx_lnd => seb__mem % latent_hflx_lnd % ptr ( ics : ice , iblk )
!$cla$w acc data create(s_air, heat_tcoef)
  ALLOCATE ( s_air ( nc ) , heat_tcoef ( nc ) )
!$acc enter data create(s_air, heat_tcoef)
!$acc data &
!$acc present(s_air(:),t_bcoef(:),sensible_hflx(:),s_star(:) &
!$acc ,sensible_hflx_lnd(:),latent_hflx_lnd(:),drag_srf(:),latent_hflx(:) &
!$acc ,t_acoef(:),heat_tcoef(:),evapopot(:),evapotrans(:),fract_snow(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_6 = 1 , size ( s_air , 1 )
   s_air ( claw_induction_6 ) = t_acoef ( claw_induction_6 ) * s_star (&
    claw_induction_6 ) + t_bcoef ( claw_induction_6 )
   heat_tcoef ( claw_induction_6 ) = heat_transfer_coef ( drag_srf (&
    claw_induction_6 ) , steplen )
   sensible_hflx ( claw_induction_6 ) = heat_tcoef ( claw_induction_6 ) * (&
    s_air ( claw_induction_6 ) - s_star ( claw_induction_6 ) )
   latent_hflx ( claw_induction_6 ) = alv * evapotrans ( claw_induction_6 ) +&
    ( als - alv ) * fract_snow ( claw_induction_6 ) * evapopot (&
    claw_induction_6 )
   sensible_hflx_lnd ( claw_induction_6 ) = sensible_hflx ( claw_induction_6 )
   latent_hflx_lnd ( claw_induction_6 ) = latent_hflx ( claw_induction_6 )
  END DO
!$acc end parallel
!$acc end data
!$cla$w acc end data
!$acc exit data delete(s_air, heat_tcoef)
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_land" //":update_surface_fluxes_land" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_surface_fluxes_land

END MODULE mo_seb_land

