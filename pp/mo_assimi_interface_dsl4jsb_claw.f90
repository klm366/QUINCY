MODULE mo_assimi_interface
 USE mo_jsb_control , ONLY: debug_on
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message
 USE mo_jsb_grid_class , ONLY: t_jsb_grid
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract , t_jsb_aggregator
 USE mo_jsb_process_class , ONLY: t_jsb_process
 USE mo_jsb_task_class , ONLY: t_jsb_process_task , t_jsb_task_options
 USE mo_jsb_process_class , ONLY: assimi_ , hydro_ , pheno_ , rad_ , a2l_ ,&
  seb_
 USE mo_assimi_config_class , ONLY: t_assimi_config
 USE mo_a2l_memory_class , ONLY: t_a2l_memory
 USE mo_assimi_memory_class , ONLY: t_assimi_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 USE mo_pheno_memory_class , ONLY: t_pheno_memory
 USE mo_rad_memory_class , ONLY: t_rad_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_assimi_interface"
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_assimilation_scaling_factors
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_assimilation_scaling_factors
  PROCEDURE , NOPASS :: aggregate => aggregate_assimilation_scaling_factors
 END TYPE tsk_assimilation_scaling_factors
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_canopy_cond_unstressed
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_canopy_cond_unstressed
  PROCEDURE , NOPASS :: aggregate => aggregate_canopy_cond_unstressed
 END TYPE tsk_canopy_cond_unstressed
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_canopy_cond_stressed
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_canopy_cond_stressed
  PROCEDURE , NOPASS :: aggregate => aggregate_canopy_cond_stressed
 END TYPE tsk_canopy_cond_stressed
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_assimilation
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_assimilation
  PROCEDURE , NOPASS :: aggregate => aggregate_assimilation
 END TYPE tsk_assimilation
 PRIVATE :: create_task_assimilation_scaling_factors
 PRIVATE :: create_task_canopy_cond_unstressed
 PRIVATE :: create_task_canopy_cond_stressed
 PRIVATE :: create_task_assimilation
 PRIVATE :: update_assimilation_scaling_factors
 PRIVATE :: aggregate_assimilation_scaling_factors
 PRIVATE :: update_canopy_cond_unstressed
 PRIVATE :: aggregate_canopy_cond_unstressed
 PRIVATE :: update_canopy_cond_stressed
 PRIVATE :: aggregate_canopy_cond_stressed
 PRIVATE :: update_assimilation
 PRIVATE :: aggregate_assimilation
 PUBLIC :: register_assimi_tasks
 INTERFACE tsk_assimilation_scaling_factors
  PROCEDURE create_task_assimilation_scaling_factors
 END INTERFACE
 INTERFACE tsk_canopy_cond_unstressed
  PROCEDURE create_task_canopy_cond_unstressed
 END INTERFACE
 INTERFACE tsk_canopy_cond_stressed
  PROCEDURE create_task_canopy_cond_stressed
 END INTERFACE
 INTERFACE tsk_assimilation
  PROCEDURE create_task_assimilation
 END INTERFACE

CONTAINS
 FUNCTION create_task_assimilation_scaling_factors ( model_id )&
  RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_assimilation_scaling_factors :: return_ptr )
  CALL return_ptr % construct ( name ="assimilation_scaling_factors" ,&
   process_id = assimi_ , owner_model_id = model_id )
 END FUNCTION create_task_assimilation_scaling_factors

 FUNCTION create_task_canopy_cond_unstressed ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_canopy_cond_unstressed :: return_ptr )
  CALL return_ptr % construct ( name ="canopy_cond_unstressed" , process_id =&
   assimi_ , owner_model_id = model_id )
 END FUNCTION create_task_canopy_cond_unstressed

 FUNCTION create_task_canopy_cond_stressed ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_canopy_cond_stressed :: return_ptr )
  CALL return_ptr % construct ( name ="canopy_cond_stressed" , process_id =&
   assimi_ , owner_model_id = model_id )
 END FUNCTION create_task_canopy_cond_stressed

 FUNCTION create_task_assimilation ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_assimilation :: return_ptr )
  CALL return_ptr % construct ( name ="assimilation" , process_id = assimi_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_assimilation

 SUBROUTINE register_assimi_tasks ( this , model_id )

  CLASS ( t_jsb_process ) , INTENT(INOUT) :: this
  INTEGER , INTENT(IN) :: model_id

  CALL this % register_task ( tsk_assimilation_scaling_factors ( model_id ) )
  CALL this % register_task ( tsk_canopy_cond_unstressed ( model_id ) )
  CALL this % register_task ( tsk_canopy_cond_stressed ( model_id ) )
  CALL this % register_task ( tsk_assimilation ( model_id ) )
 END SUBROUTINE register_assimi_tasks

 SUBROUTINE update_assimilation_scaling_factors ( tile , options )
  USE mo_psrad_orbit , ONLY: inquire_declination
  USE mo_jsb_grid , ONLY: get_grid
  USE mo_assimi_process , ONLY: calc_assimilation_scaling_factors

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_jsb_grid ) , POINTER :: grid
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: lat ( 1 : options % nc )
  REAL ( KIND= wp ) :: declination
  INTEGER :: i
  INTEGER :: ncanopy
  REAL ( KIND= wp ) , POINTER :: canopy_bound_lai ( : )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_interface" //&
  ":assimilation_scaling_factors"
  TYPE ( t_assimi_config ) , POINTER :: assimi__conf
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  REAL ( KIND= wp ) , POINTER :: assimi_scaling_fact_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( assimi_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":assimilation_scaling_factors" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( assimi_ ) % p % config )
  TYPE IS ( t_assimi_config )
   assimi__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  ncanopy = assimi__conf % ncanopy
  canopy_bound_lai ( 0 : ncanopy ) => assimi__conf % canopy_bound_lai ( 0 :&
   ncanopy )
  grid => get_grid ( model % grid_id )
  lat ( : ) = grid % lat ( ics : ice , iblk )
  lai => pheno__mem % lai % ptr ( ics : ice , iblk )
  assimi_scaling_fact_cl => assimi__mem % assimi_scaling_fact_cl % ptr ( ics :&
   ice , : , iblk )
!$acc data create(lat, canopy_bound_lai)
!$acc update device(lat, canopy_bound_lai) ! local variable
  CALL inquire_declination ( declination )
  IF ( model % lctlib ( tile % lcts ( 1 ) % lib_id ) % nitrogenscalingflag )&
   THEN
   DO i = 1 , ncanopy , 1
    CALL calc_assimilation_scaling_factors ( declination , canopy_bound_lai (&
     i - 1 ) , lai ( : ) , lat ( : ) , assimi_scaling_fact_cl ( : , i ) )
   END DO
  ELSE
!$acc data present(assimi_scaling_fact_cl(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_3 = 1 , size ( assimi_scaling_fact_cl , 1 )
    DO claw_induction_4 = 1 , size ( assimi_scaling_fact_cl , 2 )
     assimi_scaling_fact_cl ( claw_induction_3 , claw_induction_4 ) = 1._wp
    END DO
   END DO
!$acc end parallel
!$acc end data
  END IF
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":assimilation_scaling_factors" ) ,"Finished." )
  END IF
 END SUBROUTINE update_assimilation_scaling_factors

 SUBROUTINE aggregate_assimilation_scaling_factors ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_interface" //&
  ":aggregate_assimilation_scaling_factors"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":aggregate_assimilation_scaling_factors" ) ,"Starting on tile " // trim (&
    tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem %&
   assimi_scaling_fact_cl , ics , ice , iblk ,"mo_assimi_interface" //&
  ":aggregate_assimilation_scaling_factors" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":aggregate_assimilation_scaling_factors" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_assimilation_scaling_factors

 SUBROUTINE update_canopy_cond_unstressed ( tile , options )
  USE mo_assimi_process , ONLY: calc_assimilation_waterunlimited

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: icanopy
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_interface" //&
  ":update_canopy_cond_unstressed"
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_assimi_config ) , POINTER :: assimi__conf
  INTEGER :: ncanopy
  REAL ( KIND= wp ) , POINTER :: t_air ( : )
  REAL ( KIND= wp ) , POINTER :: press_srf ( : )
  REAL ( KIND= wp ) , POINTER :: par_down_mol ( : )
  REAL ( KIND= wp ) , POINTER :: apar_per_lai_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: canopy_cond_unlimited ( : )
  REAL ( KIND= wp ) , POINTER :: co2_air_mol ( : )
  REAL ( KIND= wp ) , POINTER :: assimi_scaling_fact_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: canopy_cond_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: lai_cl ( : , : )
  INTEGER :: claw_induction_1

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( assimi_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":update_canopy_cond_unstressed" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( assimi_ ) % p % config )
  TYPE IS ( t_assimi_config )
   assimi__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  t_air => a2l__mem % t_air % ptr ( ics : ice , iblk )
  press_srf => a2l__mem % press_srf % ptr ( ics : ice , iblk )
  co2_air_mol => a2l__mem % co2_air_mol % ptr ( ics : ice , iblk )
  par_down_mol => rad__mem % par_down_mol % ptr ( ics : ice , iblk )
  apar_per_lai_cl => rad__mem % apar_per_lai_cl % ptr ( ics : ice , : , iblk )
  lai_cl => rad__mem % lai_cl % ptr ( ics : ice , : , iblk )
  assimi_scaling_fact_cl => assimi__mem % assimi_scaling_fact_cl % ptr ( ics :&
   ice , : , iblk )
  canopy_cond_cl => assimi__mem % canopy_cond_cl % ptr ( ics : ice , : , iblk&
   )
  canopy_cond_unlimited => assimi__mem % canopy_cond_unlimited % ptr ( ics :&
   ice , iblk )
  ncanopy = assimi__conf % ncanopy
  DO icanopy = 1 , ncanopy , 1
   CALL calc_assimilation_waterunlimited ( model % lctlib ( tile % lcts ( 1 )&
    % lib_id ) % c4flag , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    carboxrate , model % lctlib ( tile % lcts ( 1 ) % lib_id ) % etransport ,&
    t_air ( : ) , press_srf ( : ) , par_down_mol ( : ) , apar_per_lai_cl ( : ,&
    icanopy ) , co2_air_mol ( : ) , assimi_scaling_fact_cl ( : , icanopy ) ,&
    canopy_cond_cl ( : , icanopy ) )
  END DO
!$acc data present(lai_cl(:,:),canopy_cond_cl(:,:),canopy_cond_unlimited(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_1 = 1 , size ( canopy_cond_unlimited , 1 )
   canopy_cond_unlimited ( claw_induction_1 ) = max ( 1.e-20_wp , sum (&
    canopy_cond_cl ( claw_induction_1 , : ) * lai_cl ( claw_induction_1 , : )&
    ) )
  END DO
!$acc end parallel
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":update_canopy_cond_unstressed" ) ,"Finished." )
  END IF
 END SUBROUTINE update_canopy_cond_unstressed

 SUBROUTINE aggregate_canopy_cond_unstressed ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_interface" //&
  ":aggregate_canopy_cond_unstressed"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":aggregate_canopy_cond_unstressed" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % canopy_cond_cl ,&
   ics , ice , iblk ,"mo_assimi_interface" //&
  ":aggregate_canopy_cond_unstressed" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem %&
   canopy_cond_unlimited , ics , ice , iblk ,"mo_assimi_interface" //&
  ":aggregate_canopy_cond_unstressed" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":aggregate_canopy_cond_unstressed" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_canopy_cond_unstressed

 SUBROUTINE update_canopy_cond_stressed ( tile , options )
  USE mo_hydro_process , ONLY: get_canopy_conductance
  USE mo_phy_schemes , ONLY: qsat

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_assimi_config ) , POINTER :: assimi__conf
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: canopy_cond_unlimited ( : )
  REAL ( KIND= wp ) , POINTER :: canopy_cond_limited ( : )
  REAL ( KIND= wp ) , POINTER :: canopy_cond_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: water_stress ( : )
  REAL ( KIND= wp ) , POINTER :: t ( : )
  REAL ( KIND= wp ) , POINTER :: q_air ( : )
  REAL ( KIND= wp ) , POINTER :: press_srf ( : )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: nc
  INTEGER :: icanopy
  INTEGER :: ncanopy
  LOGICAL , ALLOCATABLE :: q_air_gt_qsat_tmp ( : )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_interface" //&
  ":update_canopy_cond_stressed"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_5

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  nc = options % nc
  IF ( ( .NOT. tile % is_process_active ( assimi_ ) ) .OR. ( .NOT. tile %&
   is_vegetation ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //":update_canopy_cond_stressed"&
    ) ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( assimi_ ) % p % config )
  TYPE IS ( t_assimi_config )
   assimi__conf => cc
  END SELECT
  ncanopy = assimi__conf % ncanopy
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  q_air => a2l__mem % q_air % ptr ( ics : ice , iblk )
  press_srf => a2l__mem % press_srf % ptr ( ics : ice , iblk )
  t => seb__mem % t % ptr ( ics : ice , iblk )
  water_stress => hydro__mem % water_stress % ptr ( ics : ice , iblk )
  canopy_cond_unlimited => assimi__mem % canopy_cond_unlimited % ptr ( ics :&
   ice , iblk )
  canopy_cond_limited => assimi__mem % canopy_cond_limited % ptr ( ics : ice ,&
   iblk )
  canopy_cond_cl => assimi__mem % canopy_cond_cl % ptr ( ics : ice , : , iblk&
   )
!$cla$w acc data create(q_air_gt_qsat_tmp)
  ALLOCATE ( q_air_gt_qsat_tmp ( nc ) )
!$acc enter data create(q_air_gt_qsat_tmp)
!$acc data &
!$acc present(t(:),canopy_cond_limited(:),q_air_gt_qsat_tmp(:),q_air(:) &
!$acc ,press_srf(:),canopy_cond_unlimited(:),water_stress(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_0 = 1 , size ( q_air_gt_qsat_tmp , 1 )
   q_air_gt_qsat_tmp ( claw_induction_0 ) = q_air ( claw_induction_0 ) > qsat&
    ( t ( claw_induction_0 ) , press_srf ( claw_induction_0 ) )
   canopy_cond_limited ( claw_induction_0 ) = get_canopy_conductance (&
    canopy_cond_unlimited ( claw_induction_0 ) , water_stress (&
    claw_induction_0 ) , q_air_gt_qsat_tmp ( claw_induction_0 ) )
  END DO
!$acc end parallel
!$acc end data
!$cla$w acc end data
!$acc exit data delete(q_air_gt_qsat_tmp)
  DEALLOCATE ( q_air_gt_qsat_tmp )
  DO icanopy = 1 , ncanopy , 1
!$acc data &
!$acc present(canopy_cond_cl(:,icanopy),canopy_cond_limited(:) &
!$acc ,canopy_cond_unlimited(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_5 = 1 , size ( canopy_cond_cl , 1 )
    canopy_cond_cl ( claw_induction_5 , icanopy ) = canopy_cond_cl (&
     claw_induction_5 , icanopy ) * canopy_cond_limited ( claw_induction_5 ) /&
     canopy_cond_unlimited ( claw_induction_5 )
   END DO
!$acc end parallel
!$acc end data
  END DO
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //":update_canopy_cond_stressed"&
    ) ,"Finished." )
  END IF
 END SUBROUTINE update_canopy_cond_stressed

 SUBROUTINE aggregate_canopy_cond_stressed ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_interface" //&
  ":aggregate_canopy_cond_stressed"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":aggregate_canopy_cond_stressed" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem %&
   canopy_cond_limited , ics , ice , iblk ,"mo_assimi_interface" //&
  ":aggregate_canopy_cond_stressed" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % canopy_cond_cl ,&
   ics , ice , iblk ,"mo_assimi_interface" //":aggregate_canopy_cond_stressed"&
   )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //&
   ":aggregate_canopy_cond_stressed" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_canopy_cond_stressed

 SUBROUTINE update_assimilation ( tile , options )
  USE mo_assimi_process , ONLY: calc_assimilation_waterlimited ,&
   calc_npp_pot_rate

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: dtime
  INTEGER :: icanopy
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_interface" //&
  ":update_assimilation"
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  TYPE ( t_assimi_config ) , POINTER :: assimi__conf
  INTEGER :: ncanopy
  REAL ( KIND= wp ) , POINTER :: t_air ( : )
  REAL ( KIND= wp ) , POINTER :: press_srf ( : )
  REAL ( KIND= wp ) , POINTER :: par_down_mol ( : )
  REAL ( KIND= wp ) , POINTER :: apar_per_lai_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: co2_air_mol ( : )
  REAL ( KIND= wp ) , POINTER :: veg_fract_correction ( : )
  REAL ( KIND= wp ) , POINTER :: gross_assimilation_ca ( : )
  REAL ( KIND= wp ) , POINTER :: gross_assimilation ( : )
  REAL ( KIND= wp ) , POINTER :: dark_respiration_ca ( : )
  REAL ( KIND= wp ) , POINTER :: dark_respiration_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: dark_respiration ( : )
  REAL ( KIND= wp ) , POINTER :: npp_pot_rate_ca ( : )
  REAL ( KIND= wp ) , POINTER :: assimi_scaling_fact_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc_max ( : )
  REAL ( KIND= wp ) , POINTER :: canopy_cond_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: gross_assimilation_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: carbox_rate_max_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: e_transport_rate_max_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: carbox_rate_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: e_transport_rate_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: co2_conc_leaf_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: lai_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: day_npp_sum ( : )
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_6

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  dtime = options % dtime
  IF ( ( .NOT. tile % is_process_active ( assimi_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //":update_assimilation" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( assimi_ ) % p % config )
  TYPE IS ( t_assimi_config )
   assimi__conf => cc
  END SELECT
  ncanopy = assimi__conf % ncanopy
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  t_air => a2l__mem % t_air % ptr ( ics : ice , iblk )
  co2_air_mol => a2l__mem % co2_air_mol % ptr ( ics : ice , iblk )
  press_srf => a2l__mem % press_srf % ptr ( ics : ice , iblk )
  par_down_mol => rad__mem % par_down_mol % ptr ( ics : ice , iblk )
  apar_per_lai_cl => rad__mem % apar_per_lai_cl % ptr ( ics : ice , : , iblk )
  lai_cl => rad__mem % lai_cl % ptr ( ics : ice , : , iblk )
  veg_fract_correction => pheno__mem % veg_fract_correction % ptr ( ics : ice&
   , iblk )
  assimi_scaling_fact_cl => assimi__mem % assimi_scaling_fact_cl % ptr ( ics :&
   ice , : , iblk )
  fract_fpc_max => pheno__mem % fract_fpc_max % ptr ( ics : ice , iblk )
  canopy_cond_cl => assimi__mem % canopy_cond_cl % ptr ( ics : ice , : , iblk&
   )
  gross_assimilation_ca => assimi__mem % gross_assimilation_ca % ptr ( ics :&
   ice , iblk )
  gross_assimilation => assimi__mem % gross_assimilation % ptr ( ics : ice ,&
   iblk )
  dark_respiration_ca => assimi__mem % dark_respiration_ca % ptr ( ics : ice ,&
   iblk )
  dark_respiration => assimi__mem % dark_respiration % ptr ( ics : ice , iblk&
   )
  npp_pot_rate_ca => assimi__mem % npp_pot_rate_ca % ptr ( ics : ice , iblk )
  gross_assimilation_cl => assimi__mem % gross_assimilation_cl % ptr ( ics :&
   ice , : , iblk )
  dark_respiration_cl => assimi__mem % dark_respiration_cl % ptr ( ics : ice ,&
   : , iblk )
  carbox_rate_max_cl => assimi__mem % carbox_rate_max_cl % ptr ( ics : ice , :&
   , iblk )
  e_transport_rate_max_cl => assimi__mem % e_transport_rate_max_cl % ptr ( ics&
   : ice , : , iblk )
  e_transport_rate_cl => assimi__mem % e_transport_rate_cl % ptr ( ics : ice ,&
   : , iblk )
  carbox_rate_cl => assimi__mem % carbox_rate_cl % ptr ( ics : ice , : , iblk&
   )
  co2_conc_leaf_cl => assimi__mem % co2_conc_leaf_cl % ptr ( ics : ice , : ,&
   iblk )
  day_npp_sum => pheno__mem % day_npp_sum % ptr ( ics : ice , iblk )
  DO icanopy = 1 , ncanopy , 1
   CALL calc_assimilation_waterlimited ( model % lctlib ( tile % lcts ( 1 ) %&
    lib_id ) % c4flag , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
    carboxrate , model % lctlib ( tile % lcts ( 1 ) % lib_id ) % etransport ,&
    t_air ( : ) , press_srf ( : ) , par_down_mol ( : ) , apar_per_lai_cl ( : ,&
    icanopy ) , co2_air_mol ( : ) , assimi_scaling_fact_cl ( : , icanopy ) ,&
    co2_conc_leaf_cl ( : , icanopy ) , canopy_cond_cl ( : , icanopy ) ,&
    gross_assimilation_cl ( : , icanopy ) , dark_respiration_cl ( : , icanopy&
    ) , carbox_rate_max_cl ( : , icanopy ) , e_transport_rate_max_cl ( : ,&
    icanopy ) , carbox_rate_cl ( : , icanopy ) , e_transport_rate_cl ( : ,&
    icanopy ) )
  END DO
!$acc data &
!$acc present(dark_respiration_ca(:),dark_respiration_cl(:,:) &
!$acc ,gross_assimilation_cl(:,:),lai_cl(:,:),veg_fract_correction(:) &
!$acc ,dark_respiration(:),gross_assimilation_ca(:),fract_fpc_max(:) &
!$acc ,gross_assimilation(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_2 = 1 , size ( gross_assimilation_ca , 1 )
   gross_assimilation_ca ( claw_induction_2 ) = sum ( gross_assimilation_cl (&
    claw_induction_2 , : ) * lai_cl ( claw_induction_2 , : ) )
   gross_assimilation ( claw_induction_2 ) = gross_assimilation_ca (&
    claw_induction_2 ) * veg_fract_correction ( claw_induction_2 ) *&
    fract_fpc_max ( claw_induction_2 )
   dark_respiration_ca ( claw_induction_2 ) = sum ( dark_respiration_cl (&
    claw_induction_2 , : ) * lai_cl ( claw_induction_2 , : ) )
   dark_respiration ( claw_induction_2 ) = dark_respiration_ca (&
    claw_induction_2 ) * veg_fract_correction ( claw_induction_2 ) *&
    fract_fpc_max ( claw_induction_2 )
  END DO
!$acc end parallel
!$acc end data
  CALL calc_npp_pot_rate ( gross_assimilation_ca ( : ) , dark_respiration_ca (&
   : ) , npp_pot_rate_ca ( : ) )
!$acc data present(day_npp_sum(:),npp_pot_rate_ca(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_6 = 1 , size ( day_npp_sum , 1 )
   day_npp_sum ( claw_induction_6 ) = day_npp_sum ( claw_induction_6 ) +&
    npp_pot_rate_ca ( claw_induction_6 ) * dtime
  END DO
!$acc end parallel
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //":update_assimilation" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_assimilation

 SUBROUTINE aggregate_assimilation ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_interface" //&
  ":aggregate_assimilation"

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //":aggregate_assimilation" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
  TYPE IS ( t_assimi_memory )
   assimi__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % canopy_cond_cl ,&
   ics , ice , iblk ,"mo_assimi_interface" //":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem %&
   gross_assimilation_ca , ics , ice , iblk ,"mo_assimi_interface" //&
  ":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % gross_assimilation&
   , ics , ice , iblk ,"mo_assimi_interface" //":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem %&
   dark_respiration_ca , ics , ice , iblk ,"mo_assimi_interface" //&
  ":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % dark_respiration ,&
   ics , ice , iblk ,"mo_assimi_interface" //":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % npp_pot_rate_ca ,&
   ics , ice , iblk ,"mo_assimi_interface" //":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % co2_conc_leaf_cl ,&
   ics , ice , iblk ,"mo_assimi_interface" //":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % carbox_rate_cl ,&
   ics , ice , iblk ,"mo_assimi_interface" //":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem % carbox_rate_max_cl&
   , ics , ice , iblk ,"mo_assimi_interface" //":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem %&
   e_transport_rate_cl , ics , ice , iblk ,"mo_assimi_interface" //&
  ":aggregate_assimilation" )
  CALL weighted_by_fract % aggregate ( tile , assimi__mem %&
   e_transport_rate_max_cl , ics , ice , iblk ,"mo_assimi_interface" //&
  ":aggregate_assimilation" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_assimi_interface" //":aggregate_assimilation" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_assimilation

END MODULE mo_assimi_interface

