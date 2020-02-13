MODULE mo_turb_interface
 USE mo_jsb_control , ONLY: debug_on
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message , finish
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract , t_jsb_aggregator
 USE mo_jsb_process_class , ONLY: t_jsb_process
 USE mo_jsb_task_class , ONLY: t_jsb_process_task , t_jsb_task_options
 USE mo_jsb_process_class , ONLY: turb_ , seb_ , hydro_ , pheno_ , assimi_ ,&
  sse_ , a2l_
 USE mo_hydro_config_class , ONLY: t_hydro_config
 USE mo_seb_config_class , ONLY: t_seb_config
 USE mo_a2l_memory_class , ONLY: t_a2l_memory
 USE mo_turb_memory_class , ONLY: t_turb_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 USE mo_pheno_memory_class , ONLY: t_pheno_memory
 USE mo_assimi_memory_class , ONLY: t_assimi_memory
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_turb_interface"
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_surface_roughness
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_roughness
  PROCEDURE , NOPASS :: aggregate => aggregate_roughness
 END TYPE tsk_surface_roughness
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_humidity_scaling
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_humidity_scaling
  PROCEDURE , NOPASS :: aggregate => aggregate_humidity_scaling
 END TYPE tsk_humidity_scaling
 PRIVATE :: create_task_surface_roughness
 PRIVATE :: create_task_humidity_scaling
 PRIVATE :: update_roughness
 PRIVATE :: aggregate_roughness
 PRIVATE :: update_humidity_scaling
 PRIVATE :: aggregate_humidity_scaling
 PUBLIC :: register_turb_tasks
 INTERFACE tsk_surface_roughness
  PROCEDURE create_task_surface_roughness
 END INTERFACE
 INTERFACE tsk_humidity_scaling
  PROCEDURE create_task_humidity_scaling
 END INTERFACE

CONTAINS
 FUNCTION create_task_surface_roughness ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_surface_roughness :: return_ptr )
  CALL return_ptr % construct ( name ="surface_roughness" , process_id = turb_&
   , owner_model_id = model_id )
 END FUNCTION create_task_surface_roughness

 FUNCTION create_task_humidity_scaling ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_humidity_scaling :: return_ptr )
  CALL return_ptr % construct ( name ="humidity_scaling" , process_id = turb_&
   , owner_model_id = model_id )
 END FUNCTION create_task_humidity_scaling

 SUBROUTINE register_turb_tasks ( this , model_id )

  CLASS ( t_jsb_process ) , INTENT(INOUT) :: this
  INTEGER , INTENT(IN) :: model_id

  CALL this % register_task ( tsk_surface_roughness ( model_id ) )
  CALL this % register_task ( tsk_humidity_scaling ( model_id ) )
 END SUBROUTINE register_turb_tasks

 SUBROUTINE update_roughness ( tile , options )
  USE mo_turb_constants , ONLY: blending_height , roughness_snow

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_config ) , POINTER :: seb__conf
  TYPE ( t_turb_memory ) , POINTER :: turb__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: rough_m ( : )
  REAL ( KIND= wp ) , POINTER :: rough_h ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_soil ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_lice ( : )
  TYPE ( t_jsb_model ) , POINTER :: model
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_turb_interface" //&
  ":update_roughness"
  INTEGER :: claw_induction_6
  INTEGER :: claw_induction_7
  INTEGER :: claw_induction_8
  INTEGER :: claw_induction_9

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( turb_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_turb_interface" //":update_roughness" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( seb_ ) % p % config )
  TYPE IS ( t_seb_config )
   seb__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( turb_ ) % p )
  TYPE IS ( t_turb_memory )
   turb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  rough_m => turb__mem % rough_m % ptr ( ics : ice , iblk )
  rough_h => turb__mem % rough_h % ptr ( ics : ice , iblk )
  IF ( tile % is_process_active ( sse_ ) ) THEN
   fract_snow_soil => hydro__mem % fract_snow_soil % ptr ( ics : ice , iblk )
!$acc data present(rough_m(:),fract_snow_soil(:),rough_h(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_6 = 1 , size ( rough_h , 1 )
    rough_h ( claw_induction_6 ) = fract_snow_soil ( claw_induction_6 ) / log&
     ( blending_height / roughness_snow ) ** ( 2 ) + ( 1._wp - fract_snow_soil&
     ( claw_induction_6 ) ) / log ( blending_height / min ( 1._wp , rough_m (&
     claw_induction_6 ) ) ) ** ( 2 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
   IF ( tile % contains_lake ) THEN
    IF ( seb__conf % l_ice_on_lakes ) THEN
     fract_snow_lice => hydro__mem % fract_snow_lice % ptr ( ics : ice , iblk&
      )
!$acc data present(rough_m(:),rough_h(:),fract_snow_lice(:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_9 = 1 , size ( rough_h , 1 )
      rough_h ( claw_induction_9 ) = fract_snow_lice ( claw_induction_9 ) /&
       log ( blending_height / roughness_snow ) ** ( 2 ) + ( 1._wp -&
       fract_snow_lice ( claw_induction_9 ) ) / log ( blending_height / min (&
       1._wp , rough_m ( claw_induction_9 ) ) ) ** ( 2 )
     END DO
!$acc end parallel
!$acc end data
    ELSE
!$acc data present(rough_m(:),rough_h(:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_8 = 1 , size ( rough_h , 1 )
      rough_h ( claw_induction_8 ) = log ( blending_height / min ( 1._wp ,&
       rough_m ( claw_induction_8 ) ) ) ** ( 2 )
     END DO
!$acc end parallel
!$acc end data
    END IF
   ELSE
    CALL finish ( trim ("mo_turb_interface" //":update_roughness" ) ,&
    "Don""t know what to do." )
   END IF
  END IF
!$acc data present(rough_h(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_7 = 1 , size ( rough_h , 1 )
   rough_h ( claw_induction_7 ) = min ( 1._wp , blending_height * exp ( ( -&
    1._wp / sqrt ( rough_h ( claw_induction_7 ) ) ) ) )
  END DO
!$acc end parallel
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_turb_interface" //":update_roughness" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_roughness

 SUBROUTINE aggregate_roughness ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_turb_memory ) , POINTER :: turb__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_turb_interface" //&
  ":aggregate_roughness"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_turb_interface" //":aggregate_roughness" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( turb_ ) % p )
  TYPE IS ( t_turb_memory )
   turb__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , turb__mem % rough_h , ics , ice&
   , iblk ,"mo_turb_interface" //":aggregate_roughness" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_turb_interface" //":aggregate_roughness" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_roughness

 SUBROUTINE update_humidity_scaling ( tile , options )
  USE mo_sse_process , ONLY: relative_humidity_soil

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_config ) , POINTER :: hydro__conf
  TYPE ( t_turb_memory ) , POINTER :: turb__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  REAL ( KIND= wp ) , POINTER :: fract_fpc ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow ( : )
  REAL ( KIND= wp ) , POINTER :: fract_water ( : )
  REAL ( KIND= wp ) , POINTER :: max_moist ( : )
  REAL ( KIND= wp ) , POINTER :: w_soil_column ( : )
  REAL ( KIND= wp ) , POINTER :: fact_q_air ( : )
  REAL ( KIND= wp ) , POINTER :: fact_qsat_srf ( : )
  REAL ( KIND= wp ) , POINTER :: fact_qsat_trans_srf ( : )
  REAL ( KIND= wp ) , POINTER :: qsat_star ( : )
  REAL ( KIND= wp ) , POINTER :: canopy_cond_limited ( : )
  REAL ( KIND= wp ) , POINTER :: pch ( : )
  REAL ( KIND= wp ) , POINTER :: wind_air ( : )
  REAL ( KIND= wp ) , POINTER :: q_air ( : )
  REAL ( KIND= wp ) , POINTER :: w_soil_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_fc_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_pwp_sl ( : , : )
  REAL ( KIND= wp ) :: fact_qsat_veg ( 1 : options % nc )
  REAL ( KIND= wp ) :: fact_qsat_trans ( 1 : options % nc )
  REAL ( KIND= wp ) :: fact_qsat_soil ( 1 : options % nc )
  REAL ( KIND= wp ) :: fact_qair_veg ( 1 : options % nc )
  REAL ( KIND= wp ) :: fact_qair_soil ( 1 : options % nc )
  REAL ( KIND= wp ) :: rel_hum ( 1 : options % nc )
  REAL ( KIND= wp ) :: zfract_fpc ( 1 : options % nc )
  REAL ( KIND= wp ) :: w_soil_wilt_fract_tmp
  LOGICAL :: mask ( 1 : options % nc )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: i
  TYPE ( t_jsb_model ) , POINTER :: model
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_turb_interface" //&
  ":update_humidity_scaling"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4
  INTEGER :: claw_induction_5
  INTEGER :: claw_induction_10
  INTEGER :: claw_induction_11

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( turb_ ) ) .OR. tile % contains_lake&
   ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_turb_interface" //":update_humidity_scaling" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( hydro_ ) % p % config )
  TYPE IS ( t_hydro_config )
   hydro__conf => cc
  END SELECT
!$acc data create(fact_qsat_veg,fact_qsat_trans,fact_qsat_soil, &
!$acc fact_qair_veg,fact_qair_soil,rel_hum,zfract_fpc,mask)
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( turb_ ) % p )
  TYPE IS ( t_turb_memory )
   turb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  IF ( tile % contains_vegetation ) THEN
   SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
   TYPE IS ( t_pheno_memory )
    pheno__mem => mm
   END SELECT
   IF ( tile % is_process_active ( assimi_ ) ) THEN
    SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
    TYPE IS ( t_assimi_memory )
     assimi__mem => mm
    END SELECT
   END IF
  END IF
  pch => a2l__mem % pch % ptr ( ics : ice , iblk )
  wind_air => a2l__mem % wind_air % ptr ( ics : ice , iblk )
  q_air => a2l__mem % q_air % ptr ( ics : ice , iblk )
  fract_snow => hydro__mem % fract_snow % ptr ( ics : ice , iblk )
  fact_q_air => turb__mem % fact_q_air % ptr ( ics : ice , iblk )
  fact_qsat_srf => turb__mem % fact_qsat_srf % ptr ( ics : ice , iblk )
  fact_qsat_trans_srf => turb__mem % fact_qsat_trans_srf % ptr ( ics : ice ,&
   iblk )
  IF ( tile % contains_soil ) THEN
   qsat_star => seb__mem % qsat_star % ptr ( ics : ice , iblk )
   fract_water => hydro__mem % fract_water % ptr ( ics : ice , iblk )
   max_moist => hydro__mem % max_moist % ptr ( ics : ice , iblk )
   w_soil_column => hydro__mem % w_soil_column % ptr ( ics : ice , iblk )
   w_soil_sl => hydro__mem % w_soil_sl % ptr ( ics : ice , : , iblk )
   w_soil_fc_sl => hydro__mem % w_soil_fc_sl % ptr ( ics : ice , : , iblk )
   w_soil_pwp_sl => hydro__mem % w_soil_pwp_sl % ptr ( ics : ice , : , iblk )
  END IF
  IF ( tile % contains_vegetation ) THEN
   fract_fpc => pheno__mem % fract_fpc % ptr ( ics : ice , iblk )
   IF ( tile % is_process_active ( assimi_ ) ) THEN
    canopy_cond_limited => assimi__mem % canopy_cond_limited % ptr ( ics : ice&
     , iblk )
   ELSE
    canopy_cond_limited => hydro__mem % canopy_cond_limited % ptr ( ics : ice&
     , iblk )
   END IF
   IF ( model % config % l_compat401 ) THEN
!$acc data present(w_soil_sl(:,1),mask(:),w_soil_pwp_sl(:,1))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_11 = 1 , size ( mask , 1 )
     mask ( claw_induction_11 ) = w_soil_sl ( claw_induction_11 , 1 ) >&
      w_soil_pwp_sl ( claw_induction_11 , 1 )
    END DO
!$acc end parallel
!$acc end data
   ELSE
    w_soil_wilt_fract_tmp = hydro__conf % w_soil_wilt_fract
!$acc data present(w_soil_column(:),mask(:),max_moist(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_10 = 1 , size ( mask , 1 )
     mask ( claw_induction_10 ) = w_soil_column ( claw_induction_10 ) >&
      w_soil_wilt_fract_tmp * max_moist ( claw_induction_10 )
    END DO
!$acc end parallel
!$acc end data
   END IF
!$acc data &
!$acc present(fract_snow,fract_water,pch,wind_air,canopy_cond_limited &
!$acc ,fact_qsat_veg,fact_qsat_trans,mask)
!$acc parallel loop gang vector
   DO i = 1 , size ( mask , dim = 1 ) , 1
    IF ( mask ( i ) ) THEN
     fact_qsat_veg ( i ) = fract_snow ( i ) + ( 1._wp - fract_snow ( i ) ) * (&
      fract_water ( i ) + ( 1._wp - fract_water ( i ) ) / ( 1._wp + pch ( i )&
      * max ( 1._wp , wind_air ( i ) ) / max ( 1.e-20_wp , canopy_cond_limited&
      ( i ) ) ) )
     fact_qsat_trans ( i ) = ( 1._wp - fract_snow ( i ) ) * ( 1._wp -&
      fract_water ( i ) ) / ( 1._wp + pch ( i ) * max ( 1._wp , wind_air ( i )&
      ) / max ( 1.e-20_wp , canopy_cond_limited ( i ) ) )
    ELSE
     fact_qsat_veg ( i ) = fract_snow ( i ) + ( 1._wp - fract_snow ( i ) ) *&
      fract_water ( i )
     fact_qsat_trans ( i ) = 0._wp
    END IF
   END DO
!$acc end data
!$acc data present(fact_qsat_veg(:),fact_qair_veg(:),zfract_fpc(:),fract_fpc(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_0 = 1 , size ( fact_qair_veg , 1 )
    fact_qair_veg ( claw_induction_0 ) = fact_qsat_veg ( claw_induction_0 )
    zfract_fpc ( claw_induction_0 ) = fract_fpc ( claw_induction_0 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data &
!$acc present(fact_qsat_veg(:),fact_qair_veg(:),zfract_fpc(:) &
!$acc ,fact_qsat_trans(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_1 = 1 , size ( fact_qsat_veg , 1 )
    fact_qsat_veg ( claw_induction_1 ) = 0._wp
    fact_qsat_trans ( claw_induction_1 ) = 0._wp
    fact_qair_veg ( claw_induction_1 ) = 0._wp
    zfract_fpc ( claw_induction_1 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( tile % contains_soil ) THEN
!$acc data present(rel_hum(:),w_soil_sl(:,1),w_soil_fc_sl(:,1))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_5 = 1 , size ( rel_hum , 1 )
    rel_hum ( claw_induction_5 ) = relative_humidity_soil ( w_soil_sl (&
     claw_induction_5 , 1 ) , w_soil_fc_sl ( claw_induction_5 , 1 ) )
   END DO
!$acc end parallel
!$acc end data
!$acc data &
!$acc present(fact_qsat_soil,fact_qair_soil,qsat_star,rel_hum,q_air &
!$acc ,fract_snow,fract_water)
!$acc parallel loop gang vector
   DO i = 1 , size ( qsat_star , dim = 1 ) , 1
    IF ( qsat_star ( i ) * rel_hum ( i ) > q_air ( i ) .AND. rel_hum ( i ) >&
     1.e-10_wp ) THEN
     fact_qsat_soil ( i ) = fract_snow ( i ) + ( 1._wp - fract_snow ( i ) ) *&
      ( fract_water ( i ) + ( 1._wp - fract_water ( i ) ) * rel_hum ( i ) )
     fact_qair_soil ( i ) = 1._wp
    ELSE
     fact_qsat_soil ( i ) = fract_snow ( i ) + ( 1._wp - fract_snow ( i ) ) *&
      fract_water ( i )
     fact_qair_soil ( i ) = fact_qsat_soil ( i )
    END IF
    IF ( q_air ( i ) > qsat_star ( i ) ) THEN
     fact_qsat_soil ( i ) = 1._wp
     fact_qair_soil ( i ) = 1._wp
    END IF
   END DO
!$acc end data
  ELSE
   IF ( tile % contains_glacier ) THEN
!$acc data present(fact_qair_soil(:),fact_qsat_soil(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_2 = 1 , size ( fact_qsat_soil , 1 )
     fact_qsat_soil ( claw_induction_2 ) = 1._wp
     fact_qair_soil ( claw_induction_2 ) = 1._wp
    END DO
!$acc end parallel
!$acc end data
   ELSE
!$acc data present(fact_qair_soil(:),fact_qsat_soil(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_3 = 1 , size ( fact_qsat_soil , 1 )
     fact_qsat_soil ( claw_induction_3 ) = 0._wp
     fact_qair_soil ( claw_induction_3 ) = 0._wp
    END DO
!$acc end parallel
!$acc end data
   END IF
  END IF
!$acc data &
!$acc present(fact_qair_soil(:),fact_qsat_trans_srf(:),fact_qsat_veg(:) &
!$acc ,fact_q_air(:),zfract_fpc(:),fact_qair_veg(:),fact_qsat_soil(:) &
!$acc ,fact_qsat_trans(:),fact_qsat_srf(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_4 = 1 , size ( fact_qsat_srf , 1 )
   fact_qsat_srf ( claw_induction_4 ) = zfract_fpc ( claw_induction_4 ) *&
    fact_qsat_veg ( claw_induction_4 ) + ( 1._wp - zfract_fpc (&
    claw_induction_4 ) ) * fact_qsat_soil ( claw_induction_4 )
   fact_q_air ( claw_induction_4 ) = zfract_fpc ( claw_induction_4 ) *&
    fact_qair_veg ( claw_induction_4 ) + ( 1._wp - zfract_fpc (&
    claw_induction_4 ) ) * fact_qair_soil ( claw_induction_4 )
   fact_qsat_trans_srf ( claw_induction_4 ) = zfract_fpc ( claw_induction_4 )&
    * fact_qsat_trans ( claw_induction_4 )
  END DO
!$acc end parallel
!$acc end data
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_turb_interface" //":update_humidity_scaling" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_humidity_scaling

 SUBROUTINE aggregate_humidity_scaling ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_turb_memory ) , POINTER :: turb__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_turb_interface" //&
  ":aggregate_humidity_scaling"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_turb_interface" //":aggregate_humidity_scaling" )&
    ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( turb_ ) % p )
  TYPE IS ( t_turb_memory )
   turb__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , turb__mem % fact_q_air , ics ,&
   ice , iblk ,"mo_turb_interface" //":aggregate_humidity_scaling" )
  CALL weighted_by_fract % aggregate ( tile , turb__mem % fact_qsat_srf , ics&
   , ice , iblk ,"mo_turb_interface" //":aggregate_humidity_scaling" )
  CALL weighted_by_fract % aggregate ( tile , turb__mem % fact_qsat_trans_srf&
   , ics , ice , iblk ,"mo_turb_interface" //":aggregate_humidity_scaling" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_turb_interface" //":aggregate_humidity_scaling" )&
    ,"Finished." )
  END IF
 END SUBROUTINE aggregate_humidity_scaling

END MODULE mo_turb_interface

