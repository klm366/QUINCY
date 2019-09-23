MODULE mo_hydro_interface
 USE mo_jsb_control , ONLY: debug_on
 USE mo_jsb_time , ONLY: is_time_experiment_start
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message , finish
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_grid_class , ONLY: t_jsb_grid , t_jsb_vgrid
 USE mo_jsb_grid , ONLY: get_grid , get_vgrid
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract , t_jsb_aggregator
 USE mo_jsb_process_class , ONLY: t_jsb_process
 USE mo_jsb_task_class , ONLY: t_jsb_process_task , t_jsb_task_options
 USE mo_jsb_process_class , ONLY: seb_ , sse_ , turb_ , pheno_ , a2l_ , hydro_&
  , assimi_
 USE mo_sse_config_class , ONLY: t_sse_config
 USE mo_hydro_config_class , ONLY: t_hydro_config
 USE mo_a2l_memory_class , ONLY: t_a2l_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_sse_memory_class , ONLY: t_sse_memory
 USE mo_turb_memory_class , ONLY: t_turb_memory
 USE mo_pheno_memory_class , ONLY: t_pheno_memory
 USE mo_assimi_memory_class , ONLY: t_assimi_memory
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_hydro_interface"
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_surface_hydrology
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_surface_hydrology
  PROCEDURE , NOPASS :: aggregate => aggregate_surface_hydrology
 END TYPE tsk_surface_hydrology
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_soil_properties
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_soil_properties
  PROCEDURE , NOPASS :: aggregate => aggregate_soil_properties
 END TYPE tsk_soil_properties
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_soil_hydrology
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_soil_hydrology
  PROCEDURE , NOPASS :: aggregate => aggregate_soil_hydrology
 END TYPE tsk_soil_hydrology
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_evaporation
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_evaporation
  PROCEDURE , NOPASS :: aggregate => aggregate_evaporation
 END TYPE tsk_evaporation
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
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_snow_and_ice_hydrology
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_snow_and_ice_hydrology
  PROCEDURE , NOPASS :: aggregate => aggregate_snow_and_ice_hydrology
 END TYPE tsk_snow_and_ice_hydrology
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_snow_and_skin_fraction
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_snow_and_skin_fraction
  PROCEDURE , NOPASS :: aggregate => aggregate_snow_and_skin_fraction
 END TYPE tsk_snow_and_skin_fraction
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_water_balance
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_water_balance
  PROCEDURE , NOPASS :: aggregate => aggregate_water_balance
 END TYPE tsk_water_balance
 PRIVATE :: create_task_surface_hydrology
 PRIVATE :: create_task_soil_properties
 PRIVATE :: create_task_soil_hydrology
 PRIVATE :: create_task_canopy_cond_unstressed
 PRIVATE :: create_task_canopy_cond_stressed
 PRIVATE :: create_task_evaporation
 PRIVATE :: create_task_snow_and_ice_hydrology
 PRIVATE :: create_task_snow_and_skin_fraction
 PRIVATE :: create_task_water_balance
 PRIVATE :: update_surface_hydrology
 PRIVATE :: aggregate_surface_hydrology
 PRIVATE :: update_soil_properties
 PRIVATE :: aggregate_soil_properties
 PRIVATE :: update_soil_hydrology
 PRIVATE :: aggregate_soil_hydrology
 PRIVATE :: update_canopy_cond_unstressed
 PRIVATE :: aggregate_canopy_cond_unstressed
 PRIVATE :: update_canopy_cond_stressed
 PRIVATE :: aggregate_canopy_cond_stressed
 PRIVATE :: update_evaporation
 PRIVATE :: aggregate_evaporation
 PRIVATE :: update_snow_and_ice_hydrology
 PRIVATE :: aggregate_snow_and_ice_hydrology
 PRIVATE :: update_snow_and_skin_fraction
 PRIVATE :: aggregate_snow_and_skin_fraction
 PRIVATE :: update_water_balance
 PRIVATE :: aggregate_water_balance
 PUBLIC :: register_hydro_tasks
 INTERFACE tsk_surface_hydrology
  PROCEDURE create_task_surface_hydrology
 END INTERFACE
 INTERFACE tsk_soil_properties
  PROCEDURE create_task_soil_properties
 END INTERFACE
 INTERFACE tsk_soil_hydrology
  PROCEDURE create_task_soil_hydrology
 END INTERFACE
 INTERFACE tsk_evaporation
  PROCEDURE create_task_evaporation
 END INTERFACE
 INTERFACE tsk_canopy_cond_unstressed
  PROCEDURE create_task_canopy_cond_unstressed
 END INTERFACE
 INTERFACE tsk_canopy_cond_stressed
  PROCEDURE create_task_canopy_cond_stressed
 END INTERFACE
 INTERFACE tsk_snow_and_ice_hydrology
  PROCEDURE create_task_snow_and_ice_hydrology
 END INTERFACE
 INTERFACE tsk_snow_and_skin_fraction
  PROCEDURE create_task_snow_and_skin_fraction
 END INTERFACE
 INTERFACE tsk_water_balance
  PROCEDURE create_task_water_balance
 END INTERFACE

CONTAINS
 FUNCTION create_task_surface_hydrology ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_surface_hydrology :: return_ptr )
  CALL return_ptr % construct ( name ="surface_hydrology" , process_id =&
   hydro_ , owner_model_id = model_id )
 END FUNCTION create_task_surface_hydrology

 FUNCTION create_task_soil_properties ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_soil_properties :: return_ptr )
  CALL return_ptr % construct ( name ="soil_properties" , process_id = hydro_&
   , owner_model_id = model_id )
 END FUNCTION create_task_soil_properties

 FUNCTION create_task_soil_hydrology ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_soil_hydrology :: return_ptr )
  CALL return_ptr % construct ( name ="soil_hydrology" , process_id = hydro_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_soil_hydrology

 FUNCTION create_task_canopy_cond_unstressed ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_canopy_cond_unstressed :: return_ptr )
  CALL return_ptr % construct ( name ="canopy_cond_unstressed" , process_id =&
   hydro_ , owner_model_id = model_id )
 END FUNCTION create_task_canopy_cond_unstressed

 FUNCTION create_task_canopy_cond_stressed ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_canopy_cond_stressed :: return_ptr )
  CALL return_ptr % construct ( name ="canopy_cond_stressed" , process_id =&
   hydro_ , owner_model_id = model_id )
 END FUNCTION create_task_canopy_cond_stressed

 FUNCTION create_task_evaporation ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_evaporation :: return_ptr )
  CALL return_ptr % construct ( name ="evaporation" , process_id = hydro_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_evaporation

 FUNCTION create_task_snow_and_ice_hydrology ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_snow_and_ice_hydrology :: return_ptr )
  CALL return_ptr % construct ( name ="snow_ice_hydrology" , process_id =&
   hydro_ , owner_model_id = model_id )
 END FUNCTION create_task_snow_and_ice_hydrology

 FUNCTION create_task_snow_and_skin_fraction ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_snow_and_skin_fraction :: return_ptr )
  CALL return_ptr % construct ( name ="snow_and_skin_fraction" , process_id =&
   hydro_ , owner_model_id = model_id )
 END FUNCTION create_task_snow_and_skin_fraction

 FUNCTION create_task_water_balance ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_water_balance :: return_ptr )
  CALL return_ptr % construct ( name ="water_balance" , process_id = hydro_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_water_balance

 SUBROUTINE register_hydro_tasks ( process , model_id )

  CLASS ( t_jsb_process ) , INTENT(INOUT) :: process
  INTEGER , INTENT(IN) :: model_id

  CALL process % register_task ( tsk_surface_hydrology ( model_id ) )
  CALL process % register_task ( tsk_soil_properties ( model_id ) )
  CALL process % register_task ( tsk_soil_hydrology ( model_id ) )
  CALL process % register_task ( tsk_canopy_cond_unstressed ( model_id ) )
  CALL process % register_task ( tsk_canopy_cond_stressed ( model_id ) )
  CALL process % register_task ( tsk_evaporation ( model_id ) )
  CALL process % register_task ( tsk_snow_and_skin_fraction ( model_id ) )
  CALL process % register_task ( tsk_water_balance ( model_id ) )
 END SUBROUTINE register_hydro_tasks

 SUBROUTINE update_surface_hydrology ( tile , options )
  USE mo_hydro_process , ONLY: calc_surface_hydrology_land ,&
   calc_surface_hydrology_glacier
  USE mo_jsb_physical_constants , ONLY: dens_snow

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_config ) , POINTER :: hydro__conf
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  REAL ( KIND= wp ) , POINTER :: t_unfilt ( : )
  REAL ( KIND= wp ) , POINTER :: q_snocpymlt ( : )
  REAL ( KIND= wp ) , POINTER :: t_air ( : )
  REAL ( KIND= wp ) , POINTER :: wind_10m ( : )
  REAL ( KIND= wp ) , POINTER :: rain ( : )
  REAL ( KIND= wp ) , POINTER :: snow ( : )
  REAL ( KIND= wp ) , POINTER :: w_skin ( : )
  REAL ( KIND= wp ) , POINTER :: fract_water ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow_soil ( : )
  REAL ( KIND= wp ) , POINTER :: snow_soil_dens ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow ( : )
  REAL ( KIND= wp ) , POINTER :: snow_accum ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow_can ( : )
  REAL ( KIND= wp ) , POINTER :: evapotrans_soil ( : )
  REAL ( KIND= wp ) , POINTER :: evapo_skin ( : )
  REAL ( KIND= wp ) , POINTER :: evapo_deficit ( : )
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc_max ( : )
  REAL ( KIND= wp ) , POINTER :: hcap_grnd ( : )
  REAL ( KIND= wp ) , POINTER :: evapotrans ( : )
  REAL ( KIND= wp ) , POINTER :: trans ( : )
  REAL ( KIND= wp ) , POINTER :: evapopot ( : )
  REAL ( KIND= wp ) , POINTER :: evapo_snow ( : )
  REAL ( KIND= wp ) , POINTER :: snowmelt ( : )
  REAL ( KIND= wp ) , POINTER :: runoff ( : )
  REAL ( KIND= wp ) , POINTER :: drainage ( : )
  REAL ( KIND= wp ) , POINTER :: runoff_glac ( : )
  REAL ( KIND= wp ) , POINTER :: water_excess ( : )
  REAL ( KIND= wp ) , POINTER :: w_glac ( : )
  REAL ( KIND= wp ) :: w_skin_max ( 1 : options % nc )
  REAL ( KIND= wp ) :: w_skin_canopy_max ( 1 : options % nc )
  REAL ( KIND= wp ) :: trans_tmp ( 1 : options % nc )
  LOGICAL :: lstart
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: dtime
  REAL ( KIND= wp ) :: w_skin_max_config
  TYPE ( t_jsb_model ) , POINTER :: model
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_surface_hydrology"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_11
  INTEGER :: claw_induction_12
  INTEGER :: claw_induction_13
  INTEGER :: claw_induction_14
  INTEGER :: claw_induction_15
  INTEGER :: claw_induction_16
  INTEGER :: claw_induction_17
  INTEGER :: claw_induction_29
  INTEGER :: claw_induction_30
  INTEGER :: claw_induction_31
  INTEGER :: claw_induction_32
  INTEGER :: claw_induction_33

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  dtime = options % dtime
  IF ( ( .NOT. tile % is_process_active ( hydro_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_surface_hydrology" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( hydro_ ) % p % config )
  TYPE IS ( t_hydro_config )
   hydro__conf => cc
  END SELECT
  SELECT TYPE ( cc => model % processes ( sse_ ) % p % config )
  TYPE IS ( t_sse_config )
   sse__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  IF ( tile % is_lake ) THEN
   rain => a2l__mem % rain % ptr ( ics : ice , iblk )
   snow => a2l__mem % snow % ptr ( ics : ice , iblk )
   evapotrans => hydro__mem % evapotrans % ptr ( ics : ice , iblk )
   evapopot => hydro__mem % evapopot % ptr ( ics : ice , iblk )
   runoff => hydro__mem % runoff % ptr ( ics : ice , iblk )
   drainage => hydro__mem % drainage % ptr ( ics : ice , iblk )
   w_snow => hydro__mem % w_snow % ptr ( ics : ice , iblk )
   evapo_snow => hydro__mem % evapo_snow % ptr ( ics : ice , iblk )
!$acc data &
!$acc present(snow(:),drainage(:),rain(:),runoff(:),evapo_snow(:),evapopot(:) &
!$acc ,evapotrans(:),w_snow(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_0 = 1 , size ( runoff , 1 )
    runoff ( claw_induction_0 ) = rain ( claw_induction_0 ) + snow (&
     claw_induction_0 ) + evapotrans ( claw_induction_0 )
    drainage ( claw_induction_0 ) = 0._wp
    w_snow ( claw_induction_0 ) = 0._wp
    evapo_snow ( claw_induction_0 ) = evapopot ( claw_induction_0 )
   END DO
!$acc end parallel
!$acc end data
   RETURN
  END IF
!$acc data pcreate(w_skin_max, w_skin_canopy_max, trans_tmp)
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  IF ( tile % contains_vegetation ) THEN
   SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
   TYPE IS ( t_pheno_memory )
    pheno__mem => mm
   END SELECT
  END IF
  t_air => a2l__mem % t_air % ptr ( ics : ice , iblk )
  wind_10m => a2l__mem % wind_10m % ptr ( ics : ice , iblk )
  rain => a2l__mem % rain % ptr ( ics : ice , iblk )
  snow => a2l__mem % snow % ptr ( ics : ice , iblk )
  hcap_grnd => sse__mem % hcap_grnd % ptr ( ics : ice , iblk )
  w_snow => hydro__mem % w_snow % ptr ( ics : ice , iblk )
  w_snow_soil => hydro__mem % w_snow_soil % ptr ( ics : ice , iblk )
  snow_soil_dens => hydro__mem % snow_soil_dens % ptr ( ics : ice , iblk )
  fract_snow => hydro__mem % fract_snow % ptr ( ics : ice , iblk )
  evapotrans => hydro__mem % evapotrans % ptr ( ics : ice , iblk )
  evapopot => hydro__mem % evapopot % ptr ( ics : ice , iblk )
  evapo_snow => hydro__mem % evapo_snow % ptr ( ics : ice , iblk )
  q_snocpymlt => hydro__mem % q_snocpymlt % ptr ( ics : ice , iblk )
  snowmelt => hydro__mem % snowmelt % ptr ( ics : ice , iblk )
  runoff => hydro__mem % runoff % ptr ( ics : ice , iblk )
  drainage => hydro__mem % drainage % ptr ( ics : ice , iblk )
  t_unfilt => seb__mem % t_unfilt % ptr ( ics : ice , iblk )
  IF ( tile % is_glacier ) THEN
   w_glac => hydro__mem % w_glac % ptr ( ics : ice , iblk )
   runoff_glac => hydro__mem % runoff_glac % ptr ( ics : ice , iblk )
  ELSE
   water_excess => hydro__mem % water_excess % ptr ( ics : ice , iblk )
   fract_water => hydro__mem % fract_water % ptr ( ics : ice , iblk )
   w_skin => hydro__mem % w_skin % ptr ( ics : ice , iblk )
   snow_accum => hydro__mem % snow_accum % ptr ( ics : ice , iblk )
   evapotrans_soil => hydro__mem % evapotrans_soil % ptr ( ics : ice , iblk )
   evapo_skin => hydro__mem % evapo_skin % ptr ( ics : ice , iblk )
   evapo_deficit => hydro__mem % evapo_deficit % ptr ( ics : ice , iblk )
   IF ( tile % is_vegetation .OR. tile % is_land ) THEN
    w_snow_can => hydro__mem % w_snow_can % ptr ( ics : ice , iblk )
    trans => hydro__mem % trans % ptr ( ics : ice , iblk )
    lai => pheno__mem % lai % ptr ( ics : ice , iblk )
    fract_fpc_max => pheno__mem % fract_fpc_max % ptr ( ics : ice , iblk )
   END IF
  END IF
  IF ( is_time_experiment_start ( ) ) THEN
   lstart = .TRUE.
   w_snow ( : ) = 0._wp
  ELSE
   lstart = .FALSE.
  END IF
  w_skin_max_config = hydro__conf % w_skin_max
  IF ( tile % is_vegetation .OR. tile % is_land ) THEN
!$acc data present(lai(:),w_skin_canopy_max(:),fract_fpc_max(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_33 = 1 , size ( w_skin_canopy_max , 1 )
    w_skin_canopy_max ( claw_induction_33 ) = w_skin_max_config * lai (&
     claw_induction_33 ) * fract_fpc_max ( claw_induction_33 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(w_skin_canopy_max(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_32 = 1 , size ( w_skin_canopy_max , 1 )
    w_skin_canopy_max ( claw_induction_32 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( tile % contains_soil ) THEN
!$acc data present(w_skin_max(:),w_skin_canopy_max(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_31 = 1 , size ( w_skin_max , 1 )
    w_skin_max ( claw_induction_31 ) = w_skin_max_config + w_skin_canopy_max (&
     claw_induction_31 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(w_skin_max(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_30 = 1 , size ( w_skin_max , 1 )
    w_skin_max ( claw_induction_30 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( tile % is_glacier ) THEN
   CALL calc_surface_hydrology_glacier ( lstart , dtime , t_unfilt ( : ) ,&
    fract_snow ( : ) , hcap_grnd ( : ) , evapotrans ( : ) , evapopot ( : ) ,&
    rain ( : ) , snow ( : ) , w_glac ( : ) , q_snocpymlt ( : ) , snowmelt ( :&
    ) , runoff_glac ( : ) , runoff ( : ) )
!$acc data present(drainage(:),w_snow_soil(:),evapo_snow(:),evapopot(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_1 = 1 , size ( drainage , 1 )
    drainage ( claw_induction_1 ) = 0._wp
    w_snow_soil ( claw_induction_1 ) = 0._wp
    evapo_snow ( claw_induction_1 ) = evapopot ( claw_induction_1 )
   END DO
!$acc end parallel
!$acc end data
   IF ( model % config % l_compat401 ) THEN
!$acc data present(snow_soil_dens(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_17 = 1 , size ( snow_soil_dens , 1 )
     snow_soil_dens ( claw_induction_17 ) = 330._wp
    END DO
!$acc end parallel
!$acc end data
   ELSE
!$acc data present(snow_soil_dens(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_16 = 1 , size ( snow_soil_dens , 1 )
     snow_soil_dens ( claw_induction_16 ) = dens_snow
    END DO
!$acc end parallel
!$acc end data
   END IF
  ELSE
   IF ( tile % is_vegetation ) THEN
!$acc data present(trans(:),trans_tmp(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_15 = 1 , size ( trans_tmp , 1 )
     trans_tmp ( claw_induction_15 ) = trans ( claw_induction_15 )
    END DO
!$acc end parallel
!$acc end data
   ELSE
!$acc data present(trans_tmp(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_14 = 1 , size ( trans_tmp , 1 )
     trans_tmp ( claw_induction_14 ) = 0._wp
    END DO
!$acc end parallel
!$acc end data
   END IF
   CALL calc_surface_hydrology_land ( lstart , dtime , sse__conf % l_dynsnow ,&
    t_unfilt ( : ) , wind_10m ( : ) , t_air ( : ) , w_skin_canopy_max ( : ) ,&
    w_skin_max ( : ) , fract_snow ( : ) , fract_water ( : ) , hcap_grnd ( : )&
    , evapotrans ( : ) , evapopot ( : ) , trans_tmp ( : ) , rain ( : ) , snow&
    ( : ) , w_skin ( : ) , w_snow_soil ( : ) , w_snow_can ( : ) ,&
    snow_soil_dens ( : ) , q_snocpymlt ( : ) , snow_accum ( : ) , snowmelt ( :&
    ) , evapotrans_soil ( : ) , evapo_skin ( : ) , evapo_snow ( : ) ,&
    water_excess ( : ) , evapo_deficit ( : ) )
   IF ( model % config % l_compat401 ) THEN
!$acc data present(snow_soil_dens(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_13 = 1 , size ( snow_soil_dens , 1 )
     snow_soil_dens ( claw_induction_13 ) = 330._wp
    END DO
!$acc end parallel
!$acc end data
   END IF
  END IF
  IF ( tile % is_glacier ) THEN
!$acc data present(w_snow(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_29 = 1 , size ( w_snow , 1 )
    w_snow ( claw_induction_29 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  ELSE
   IF ( tile % is_bare ) THEN
!$acc data present(w_snow_soil(:),w_snow(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_12 = 1 , size ( w_snow , 1 )
     w_snow ( claw_induction_12 ) = w_snow_soil ( claw_induction_12 )
    END DO
!$acc end parallel
!$acc end data
   ELSE
    IF ( tile % is_vegetation .OR. tile % is_land ) THEN
!$acc data present(w_snow_soil(:),w_snow_can(:),w_snow(:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_11 = 1 , size ( w_snow , 1 )
      w_snow ( claw_induction_11 ) = w_snow_soil ( claw_induction_11 ) +&
       w_snow_can ( claw_induction_11 )
     END DO
!$acc end parallel
!$acc end data
    END IF
   END IF
  END IF
!$acc end data
 END SUBROUTINE update_surface_hydrology

 SUBROUTINE aggregate_surface_hydrology ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_surface_hydrology"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_surface_hydrology"&
    ) ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_skin , ics , ice&
   , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_snow , ics , ice&
   , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_snow_soil , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_snow_can , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % snowmelt , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % q_snocpymlt , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_glac , ics , ice&
   , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_snow_soil , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_snow_can , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % snow_accum , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % snow_soil_dens ,&
   ics , ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapotrans_soil ,&
   ics , ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapo_skin , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapo_snow , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % water_excess , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapo_deficit , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % runoff , ics , ice&
   , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % drainage , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_surface_hydrology" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_surface_hydrology"&
    ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_surface_hydrology

 SUBROUTINE update_soil_properties ( tile , options )
  USE mo_hydro_constants , ONLY: vol_porosity_org , hyd_cond_sat_org ,&
   bclapp_org , matrix_pot_org , pore_size_index_org , vol_field_cap_org ,&
   vol_p_wilt_org

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_config ) , POINTER :: hydro__conf
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: hyd_cond_sat ( : )
  REAL ( KIND= wp ) , POINTER :: vol_porosity ( : )
  REAL ( KIND= wp ) , POINTER :: bclapp ( : )
  REAL ( KIND= wp ) , POINTER :: matrix_pot ( : )
  REAL ( KIND= wp ) , POINTER :: pore_size_index ( : )
  REAL ( KIND= wp ) , POINTER :: vol_field_cap ( : )
  REAL ( KIND= wp ) , POINTER :: vol_p_wilt ( : )
  REAL ( KIND= wp ) , POINTER :: fract_org_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: hyd_cond_sat_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_porosity_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: bclapp_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: matrix_pot_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: pore_size_index_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_field_cap_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_p_wilt_sl ( : , : )
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_jsb_vgrid ) , POINTER :: soil_w
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: i
  INTEGER :: nsoil
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_soil_properties"
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( hydro_ ) ) ) THEN
   RETURN
  END IF
  IF ( tile % is_lake .OR. tile % is_glacier ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_soil_properties" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  soil_w => get_vgrid ("soil_depth_water" )
  nsoil = soil_w % n_levels
  SELECT TYPE ( cc => model % processes ( hydro_ ) % p % config )
  TYPE IS ( t_hydro_config )
   hydro__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  hyd_cond_sat => hydro__mem % hyd_cond_sat % ptr ( ics : ice , iblk )
  vol_porosity => hydro__mem % vol_porosity % ptr ( ics : ice , iblk )
  bclapp => hydro__mem % bclapp % ptr ( ics : ice , iblk )
  matrix_pot => hydro__mem % matrix_pot % ptr ( ics : ice , iblk )
  pore_size_index => hydro__mem % pore_size_index % ptr ( ics : ice , iblk )
  vol_field_cap => hydro__mem % vol_field_cap % ptr ( ics : ice , iblk )
  vol_p_wilt => hydro__mem % vol_p_wilt % ptr ( ics : ice , iblk )
  IF ( hydro__conf % l_organic ) THEN
   fract_org_sl => hydro__mem % fract_org_sl % ptr ( ics : ice , : , iblk )
  END IF
  hyd_cond_sat_sl => hydro__mem % hyd_cond_sat_sl % ptr ( ics : ice , : , iblk&
   )
  vol_porosity_sl => hydro__mem % vol_porosity_sl % ptr ( ics : ice , : , iblk&
   )
  bclapp_sl => hydro__mem % bclapp_sl % ptr ( ics : ice , : , iblk )
  matrix_pot_sl => hydro__mem % matrix_pot_sl % ptr ( ics : ice , : , iblk )
  pore_size_index_sl => hydro__mem % pore_size_index_sl % ptr ( ics : ice , :&
   , iblk )
  vol_field_cap_sl => hydro__mem % vol_field_cap_sl % ptr ( ics : ice , : ,&
   iblk )
  vol_p_wilt_sl => hydro__mem % vol_p_wilt_sl % ptr ( ics : ice , : , iblk )
  DO i = 1 , nsoil , 1
!$acc data &
!$acc present(pore_size_index(:),pore_size_index_sl(:,i),vol_porosity(:) &
!$acc ,vol_p_wilt(:),vol_porosity_sl(:,i),bclapp(:),vol_field_cap(:) &
!$acc ,bclapp_sl(:,i),matrix_pot(:),matrix_pot_sl(:,i),vol_p_wilt_sl(:,i) &
!$acc ,hyd_cond_sat(:),hyd_cond_sat_sl(:,i),vol_field_cap_sl(:,i))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_2 = 1 , size ( hyd_cond_sat_sl , 1 )
    hyd_cond_sat_sl ( claw_induction_2 , i ) = hyd_cond_sat ( claw_induction_2&
     )
    vol_porosity_sl ( claw_induction_2 , i ) = vol_porosity ( claw_induction_2&
     )
    bclapp_sl ( claw_induction_2 , i ) = bclapp ( claw_induction_2 )
    matrix_pot_sl ( claw_induction_2 , i ) = matrix_pot ( claw_induction_2 )
    pore_size_index_sl ( claw_induction_2 , i ) = pore_size_index (&
     claw_induction_2 )
    vol_field_cap_sl ( claw_induction_2 , i ) = vol_field_cap (&
     claw_induction_2 )
    vol_p_wilt_sl ( claw_induction_2 , i ) = vol_p_wilt ( claw_induction_2 )
   END DO
!$acc end parallel
!$acc end data
  END DO
  IF ( hydro__conf % l_organic ) THEN
!$acc data &
!$acc present(pore_size_index_sl(:,:),vol_porosity_sl(:,:),bclapp_sl(:,:) &
!$acc ,matrix_pot_sl(:,:),vol_field_cap_sl(:,:),vol_p_wilt_sl(:,:) &
!$acc ,fract_org_sl(:,:),hyd_cond_sat_sl(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_3 = 1 , size ( hyd_cond_sat_sl , 1 )
    DO claw_induction_4 = 1 , size ( hyd_cond_sat_sl , 2 )
     hyd_cond_sat_sl ( claw_induction_3 , claw_induction_4 ) = ( 1 -&
      fract_org_sl ( claw_induction_3 , claw_induction_4 ) ) * hyd_cond_sat_sl&
      ( claw_induction_3 , claw_induction_4 ) + fract_org_sl (&
      claw_induction_3 , claw_induction_4 ) * hyd_cond_sat_org
     vol_porosity_sl ( claw_induction_3 , claw_induction_4 ) = ( 1 -&
      fract_org_sl ( claw_induction_3 , claw_induction_4 ) ) * vol_porosity_sl&
      ( claw_induction_3 , claw_induction_4 ) + fract_org_sl (&
      claw_induction_3 , claw_induction_4 ) * vol_porosity_org
     bclapp_sl ( claw_induction_3 , claw_induction_4 ) = ( 1 - fract_org_sl (&
      claw_induction_3 , claw_induction_4 ) ) * bclapp_sl ( claw_induction_3 ,&
      claw_induction_4 ) + fract_org_sl ( claw_induction_3 , claw_induction_4&
      ) * bclapp_org
     matrix_pot_sl ( claw_induction_3 , claw_induction_4 ) = ( 1 -&
      fract_org_sl ( claw_induction_3 , claw_induction_4 ) ) * matrix_pot_sl (&
      claw_induction_3 , claw_induction_4 ) + fract_org_sl ( claw_induction_3&
      , claw_induction_4 ) * matrix_pot_org
     pore_size_index_sl ( claw_induction_3 , claw_induction_4 ) = ( 1 -&
      fract_org_sl ( claw_induction_3 , claw_induction_4 ) ) *&
      pore_size_index_sl ( claw_induction_3 , claw_induction_4 ) +&
      fract_org_sl ( claw_induction_3 , claw_induction_4 ) *&
      pore_size_index_org
     vol_field_cap_sl ( claw_induction_3 , claw_induction_4 ) = ( 1 -&
      fract_org_sl ( claw_induction_3 , claw_induction_4 ) ) *&
      vol_field_cap_sl ( claw_induction_3 , claw_induction_4 ) + fract_org_sl&
      ( claw_induction_3 , claw_induction_4 ) * vol_field_cap_org
     vol_p_wilt_sl ( claw_induction_3 , claw_induction_4 ) = ( 1 -&
      fract_org_sl ( claw_induction_3 , claw_induction_4 ) ) * vol_p_wilt_sl (&
      claw_induction_3 , claw_induction_4 ) + fract_org_sl ( claw_induction_3&
      , claw_induction_4 ) * vol_p_wilt_org
    END DO
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_soil_properties" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_soil_properties

 SUBROUTINE aggregate_soil_properties ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_soil_properties"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_soil_properties" )&
    ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % hyd_cond_sat_sl ,&
   ics , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_properties" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % vol_porosity_sl ,&
   ics , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_properties" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % bclapp_sl , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_soil_properties" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % matrix_pot_sl , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_properties" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % pore_size_index_sl&
   , ics , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_properties" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % vol_field_cap_sl ,&
   ics , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_properties" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % vol_p_wilt_sl , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_properties" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_soil_properties" )&
    ,"Finished." )
  END IF
 END SUBROUTINE aggregate_soil_properties

 SUBROUTINE update_soil_hydrology ( tile , options )
  USE mo_hydro_process , ONLY: calc_soil_hydrology

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_config ) , POINTER :: hydro__conf
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: fract_org_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: t_soil ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_rel ( : )
  REAL ( KIND= wp ) , POINTER :: w_soil_column ( : )
  REAL ( KIND= wp ) , POINTER :: max_moist ( : )
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_jsb_grid ) , POINTER :: grid
  TYPE ( t_jsb_vgrid ) , POINTER :: soil_w
  LOGICAL :: ltpe_open
  LOGICAL :: ltpe_closed
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: nc
  INTEGER :: i
  REAL ( KIND= wp ) :: dtime
  INTEGER :: nsoil
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_soil_hydrology"
  INTEGER :: claw_induction_27
  INTEGER :: claw_induction_28

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  nc = options % nc
  dtime = options % dtime
  IF ( ( .NOT. tile % is_process_active ( hydro_ ) ) ) THEN
   RETURN
  END IF
  IF ( tile % is_lake ) THEN
   RETURN
  END IF
  IF ( tile % is_glacier ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_soil_hydrology" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  grid => get_grid ( model % grid_id )
  soil_w => get_vgrid ("soil_depth_water" )
  nsoil = soil_w % n_levels
  SELECT TYPE ( cc => model % processes ( hydro_ ) % p % config )
  TYPE IS ( t_hydro_config )
   hydro__conf => cc
  END SELECT
  ltpe_open = .FALSE.
  ltpe_closed = .FALSE.
  IF ( model % config % tpe_scheme =="open" ) THEN
   ltpe_open = .TRUE.
  END IF
  IF ( model % config % tpe_scheme =="closed" ) THEN
   ltpe_closed = .TRUE.
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  w_soil_rel => hydro__mem % w_soil_rel % ptr ( ics : ice , iblk )
  w_soil_column => hydro__mem % w_soil_column % ptr ( ics : ice , iblk )
  max_moist => hydro__mem % max_moist % ptr ( ics : ice , iblk )
  t_soil => sse__mem % t_soil % ptr ( ics : ice , : , iblk )
  IF ( hydro__conf % l_organic ) THEN
   fract_org_sl => hydro__mem % fract_org_sl % ptr ( ics : ice , : , iblk )
  ELSE
   ALLOCATE ( fract_org_sl ( nc , nsoil ) )
!$acc enter data create(fract_org_sl)
!$acc data present(fract_org_sl(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_27 = 1 , size ( fract_org_sl , 1 )
    DO claw_induction_28 = 1 , size ( fract_org_sl , 2 )
     fract_org_sl ( claw_induction_27 , claw_induction_28 ) = 0._wp
    END DO
   END DO
!$acc end parallel
!$acc end data
  END IF
!$acc update device(soil_w%dz) ! todo move if possible
  CALL calc_soil_hydrology ( nc , nsoil , soil_w % dz ( : ) , grid % nlat_g ,&
   dtime , ltpe_closed , ltpe_open , hydro__mem % soil_depth_sl % ptr ( ics :&
   ice , : , iblk ) , hydro__mem % root_depth_sl % ptr ( ics : ice , : , iblk&
   ) , fract_org_sl ( : , : ) , hydro__mem % hyd_cond_sat_sl % ptr ( ics : ice&
   , : , iblk ) , hydro__mem % matrix_pot_sl % ptr ( ics : ice , : , iblk ) ,&
   hydro__mem % bclapp_sl % ptr ( ics : ice , : , iblk ) , hydro__mem %&
   pore_size_index_sl % ptr ( ics : ice , : , iblk ) , hydro__mem %&
   vol_porosity_sl % ptr ( ics : ice , : , iblk ) , hydro__mem %&
   vol_field_cap_sl % ptr ( ics : ice , : , iblk ) , hydro__mem %&
   vol_p_wilt_sl % ptr ( ics : ice , : , iblk ) , hydro__mem % max_moist % ptr&
   ( ics : ice , iblk ) , hydro__mem % oro_stddev % ptr ( ics : ice , iblk ) ,&
   t_soil ( : , 1 ) , hydro__mem % water_excess % ptr ( ics : ice , iblk ) ,&
   hydro__mem % evapotrans_soil % ptr ( ics : ice , iblk ) , hydro__mem %&
   trans % ptr ( ics : ice , iblk ) , hydro__mem % w_soil_column % ptr ( ics :&
   ice , iblk ) , hydro__mem % w_soil_sl % ptr ( ics : ice , : , iblk ) ,&
   hydro__mem % w_ice_sl % ptr ( ics : ice , : , iblk ) , hydro__mem %&
   w_soil_overflow % ptr ( ics : ice , iblk ) , hydro__mem % w_soil_sat_sl %&
   ptr ( ics : ice , : , iblk ) , hydro__mem % w_soil_fc_sl % ptr ( ics : ice&
   , : , iblk ) , hydro__mem % w_soil_pwp_sl % ptr ( ics : ice , : , iblk ) ,&
   hydro__mem % runoff % ptr ( ics : ice , iblk ) , hydro__mem % drainage %&
   ptr ( ics : ice , iblk ) )
!$acc data present(max_moist, w_soil_rel, w_soil_column)
!$acc parallel loop gang vector
  DO i = ics , ice , 1
   IF ( max_moist ( i ) > 0._wp ) THEN
    w_soil_rel ( i ) = min ( w_soil_column ( i ) / max_moist ( i ) , 1._wp )
   ELSE
    w_soil_rel ( i ) = 0._wp
   END IF
  END DO
!$acc end data
  IF ( ( .NOT. hydro__conf % l_organic ) ) THEN
!$acc exit data delete(fract_org_sl)
   DEALLOCATE ( fract_org_sl )
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_soil_hydrology" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_soil_hydrology

 SUBROUTINE aggregate_soil_hydrology ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_soil_hydrology"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_soil_hydrology" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_column , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_rel , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_sl , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_ice_sl , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_sat_sl , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_fc_sl , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_pwp_sl , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % runoff , ics , ice&
   , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % drainage , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_soil_hydrology" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_soil_hydrology" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_soil_hydrology

 SUBROUTINE update_canopy_cond_unstressed ( tile , options )
  USE mo_hydro_process , ONLY: get_canopy_conductance

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: swpar_srf_down ( : )
  REAL ( KIND= wp ) , POINTER :: canopy_cond ( : )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_canopy_cond_unstressed"
  INTEGER :: claw_induction_36

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( hydro_ ) ) .OR. ( .NOT. tile %&
   is_vegetation ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":update_canopy_cond_unstressed" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  swpar_srf_down => a2l__mem % swpar_srf_down % ptr ( ics : ice , iblk )
  lai => pheno__mem % lai % ptr ( ics : ice , iblk )
  canopy_cond => hydro__mem % canopy_cond % ptr ( ics : ice , iblk )
!$acc data present(lai(:),swpar_srf_down(:),canopy_cond(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_36 = 1 , size ( canopy_cond , 1 )
   canopy_cond ( claw_induction_36 ) = get_canopy_conductance ( lai (&
    claw_induction_36 ) , swpar_srf_down ( claw_induction_36 ) )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE update_canopy_cond_unstressed

 SUBROUTINE aggregate_canopy_cond_unstressed ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_canopy_cond_unstressed"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":aggregate_canopy_cond_unstressed" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % canopy_cond , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_canopy_cond_unstressed" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":aggregate_canopy_cond_unstressed" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_canopy_cond_unstressed

 SUBROUTINE update_canopy_cond_stressed ( tile , options )
  USE mo_hydro_process , ONLY: get_canopy_conductance ,&
   get_water_stress_factor
  USE mo_hydro_util , ONLY: get_water_in_root_zone
  USE mo_phy_schemes , ONLY: qsat

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_hydro_config ) , POINTER :: hydro__conf
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_assimi_memory ) , POINTER :: assimi__mem
  REAL ( KIND= wp ) , POINTER :: canopy_cond ( : )
  REAL ( KIND= wp ) , POINTER :: canopy_cond_limited ( : )
  REAL ( KIND= wp ) , POINTER :: soil_depth_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: root_depth_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_column ( : )
  REAL ( KIND= wp ) , POINTER :: w_soil_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_fc_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_pwp_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_root ( : )
  REAL ( KIND= wp ) , POINTER :: w_soil_root_fc ( : )
  REAL ( KIND= wp ) , POINTER :: w_soil_root_pwp ( : )
  REAL ( KIND= wp ) , POINTER :: max_moist ( : )
  REAL ( KIND= wp ) , POINTER :: water_stress ( : )
  REAL ( KIND= wp ) , POINTER :: t ( : )
  REAL ( KIND= wp ) , POINTER :: q_air ( : )
  REAL ( KIND= wp ) , POINTER :: press_srf ( : )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: config_w_soil_crit_fract
  REAL ( KIND= wp ) :: config_w_soil_wilt_fract
  LOGICAL :: q_air_gt_qsat_tmp ( 1 : options % nc )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_canopy_cond_stressed"
  INTEGER :: claw_induction_5

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( hydro_ ) ) .OR. ( .NOT. tile %&
   is_vegetation ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_canopy_cond_stressed"&
    ) ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( hydro_ ) % p % config )
  TYPE IS ( t_hydro_config )
   hydro__conf => cc
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
  IF ( tile % is_process_active ( assimi_ ) ) THEN
   SELECT TYPE ( mm => tile % mem ( assimi_ ) % p )
   TYPE IS ( t_assimi_memory )
    assimi__mem => mm
   END SELECT
  END IF
  q_air => a2l__mem % q_air % ptr ( ics : ice , iblk )
  press_srf => a2l__mem % press_srf % ptr ( ics : ice , iblk )
  t => seb__mem % t % ptr ( ics : ice , iblk )
  soil_depth_sl => hydro__mem % soil_depth_sl % ptr ( ics : ice , : , iblk )
  w_soil_sl => hydro__mem % w_soil_sl % ptr ( ics : ice , : , iblk )
  w_soil_fc_sl => hydro__mem % w_soil_fc_sl % ptr ( ics : ice , : , iblk )
  w_soil_pwp_sl => hydro__mem % w_soil_pwp_sl % ptr ( ics : ice , : , iblk )
  max_moist => hydro__mem % max_moist % ptr ( ics : ice , iblk )
  w_soil_column => hydro__mem % w_soil_column % ptr ( ics : ice , iblk )
  root_depth_sl => hydro__mem % root_depth_sl % ptr ( ics : ice , : , iblk )
  IF ( tile % is_process_active ( assimi_ ) ) THEN
   canopy_cond => assimi__mem % canopy_cond % ptr ( ics : ice , iblk )
  ELSE
   canopy_cond => hydro__mem % canopy_cond % ptr ( ics : ice , iblk )
  END IF
  canopy_cond_limited => hydro__mem % canopy_cond_limited % ptr ( ics : ice ,&
   iblk )
  w_soil_root => hydro__mem % w_soil_root % ptr ( ics : ice , iblk )
  w_soil_root_fc => hydro__mem % w_soil_root_fc % ptr ( ics : ice , iblk )
  w_soil_root_pwp => hydro__mem % w_soil_root_pwp % ptr ( ics : ice , iblk )
  water_stress => hydro__mem % water_stress % ptr ( ics : ice , iblk )
  config_w_soil_crit_fract = hydro__conf % w_soil_crit_fract
  config_w_soil_wilt_fract = hydro__conf % w_soil_wilt_fract
!$acc data present(w_soil_root, w_soil_sl, soil_depth_sl, root_depth_sl, &
!$acc w_soil_root_fc, w_soil_fc_sl, w_soil_root_pwp, w_soil_pwp_sl, w_soil_pwp_sl)
!$acc update host(w_soil_sl, soil_depth_sl, root_depth_sl, w_soil_fc_sl)
  w_soil_root ( : ) = get_water_in_root_zone ( w_soil_sl ( : , : ) ,&
   soil_depth_sl ( : , : ) , root_depth_sl ( : , : ) )
  w_soil_root_fc ( : ) = get_water_in_root_zone ( w_soil_fc_sl ( : , : ) ,&
   soil_depth_sl ( : , : ) , root_depth_sl ( : , : ) )
  w_soil_root_pwp ( : ) = get_water_in_root_zone ( w_soil_pwp_sl ( : , : ) ,&
   soil_depth_sl ( : , : ) , root_depth_sl ( : , : ) )
!$acc update device(w_soil_root,w_soil_root_fc,w_soil_root_pwp) ! computed on cpu
!$acc end data
!$acc data create(q_air_gt_qsat_tmp)
!$acc data &
!$acc present(w_soil_column(:),t(:),canopy_cond_limited(:) &
!$acc ,q_air_gt_qsat_tmp(:),max_moist(:),q_air(:),press_srf(:) &
!$acc ,water_stress(:),canopy_cond(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_5 = 1 , size ( water_stress , 1 )
   water_stress ( claw_induction_5 ) = get_water_stress_factor ( w_soil_column&
    ( claw_induction_5 ) , max_moist ( claw_induction_5 ) ,&
    config_w_soil_crit_fract , config_w_soil_wilt_fract )
   q_air_gt_qsat_tmp ( claw_induction_5 ) = q_air ( claw_induction_5 ) > qsat&
    ( t ( claw_induction_5 ) , press_srf ( claw_induction_5 ) )
   canopy_cond_limited ( claw_induction_5 ) = get_canopy_conductance (&
    canopy_cond ( claw_induction_5 ) , water_stress ( claw_induction_5 ) ,&
    q_air_gt_qsat_tmp ( claw_induction_5 ) )
  END DO
!$acc end parallel
!$acc end data
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_canopy_cond_stressed"&
    ) ,"Finished." )
  END IF
 END SUBROUTINE update_canopy_cond_stressed

 SUBROUTINE aggregate_canopy_cond_stressed ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_canopy_cond_stressed"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":aggregate_canopy_cond_stressed" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % canopy_cond_limited&
   , ics , ice , iblk ,"mo_hydro_interface" //&
  ":aggregate_canopy_cond_stressed" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_root , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_canopy_cond_stressed" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_root_fc ,&
   ics , ice , iblk ,"mo_hydro_interface" //":aggregate_canopy_cond_stressed"&
   )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_root_pwp ,&
   ics , ice , iblk ,"mo_hydro_interface" //":aggregate_canopy_cond_stressed"&
   )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % water_stress , ics&
   , ice , iblk ,"mo_hydro_interface" //":aggregate_canopy_cond_stressed" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":aggregate_canopy_cond_stressed" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_canopy_cond_stressed

 SUBROUTINE update_evaporation ( tile , options )
  USE mo_phy_schemes , ONLY: q_effective , heat_transfer_coef

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_turb_memory ) , POINTER :: turb__mem
  REAL ( KIND= wp ) , POINTER :: q_acoef ( : )
  REAL ( KIND= wp ) , POINTER :: q_bcoef ( : )
  REAL ( KIND= wp ) , POINTER :: q_acoef_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: q_bcoef_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: q_acoef_ice ( : )
  REAL ( KIND= wp ) , POINTER :: q_bcoef_ice ( : )
  REAL ( KIND= wp ) , POINTER :: drag_srf ( : )
  REAL ( KIND= wp ) , POINTER :: drag_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: drag_ice ( : )
  REAL ( KIND= wp ) , POINTER :: qsat_star ( : )
  REAL ( KIND= wp ) , POINTER :: qsat_lwtr ( : )
  REAL ( KIND= wp ) , POINTER :: qsat_lice ( : )
  REAL ( KIND= wp ) , POINTER :: fact_qsat_srf ( : )
  REAL ( KIND= wp ) , POINTER :: fact_qsat_trans_srf ( : )
  REAL ( KIND= wp ) , POINTER :: fact_q_air ( : )
  REAL ( KIND= wp ) , POINTER :: fract_lice ( : )
  REAL ( KIND= wp ) , POINTER :: evapotrans_lnd ( : )
  REAL ( KIND= wp ) , POINTER :: evapo_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: evapo_ice ( : )
  REAL ( KIND= wp ) , POINTER :: evapopot ( : )
  REAL ( KIND= wp ) , POINTER :: evapotrans ( : )
  REAL ( KIND= wp ) , POINTER :: trans ( : )
  REAL ( KIND= wp ) :: q_air ( 1 : options % nc )
  REAL ( KIND= wp ) :: q_air_eff ( 1 : options % nc )
  REAL ( KIND= wp ) :: qsat_srf_eff ( 1 : options % nc )
  REAL ( KIND= wp ) :: heat_tcoef ( 1 : options % nc )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: nc
  REAL ( KIND= wp ) :: steplen
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_evaporation"
  INTEGER :: claw_induction_6
  INTEGER :: claw_induction_7
  INTEGER :: claw_induction_10

  IF ( ( .NOT. tile % is_process_active ( hydro_ ) ) ) THEN
   RETURN
  END IF
  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  nc = options % nc
  steplen = options % steplen
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_evaporation" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
!$acc data create(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( turb_ ) % p )
  TYPE IS ( t_turb_memory )
   turb__mem => mm
  END SELECT
  evapopot => hydro__mem % evapopot % ptr ( ics : ice , iblk )
  evapotrans => hydro__mem % evapotrans % ptr ( ics : ice , iblk )
  IF ( tile % contains_lake ) THEN
   q_acoef_wtr => a2l__mem % q_acoef_wtr % ptr ( ics : ice , iblk )
   q_bcoef_wtr => a2l__mem % q_bcoef_wtr % ptr ( ics : ice , iblk )
   q_acoef_ice => a2l__mem % q_acoef_ice % ptr ( ics : ice , iblk )
   q_bcoef_ice => a2l__mem % q_bcoef_ice % ptr ( ics : ice , iblk )
   drag_wtr => a2l__mem % drag_wtr % ptr ( ics : ice , iblk )
   drag_ice => a2l__mem % drag_ice % ptr ( ics : ice , iblk )
   qsat_lwtr => seb__mem % qsat_lwtr % ptr ( ics : ice , iblk )
   qsat_lice => seb__mem % qsat_lice % ptr ( ics : ice , iblk )
   fract_lice => seb__mem % fract_lice % ptr ( ics : ice , iblk )
   evapo_wtr => hydro__mem % evapo_wtr % ptr ( ics : ice , iblk )
   evapo_ice => hydro__mem % evapo_ice % ptr ( ics : ice , iblk )
   evapopot => hydro__mem % evapopot % ptr ( ics : ice , iblk )
   evapotrans => hydro__mem % evapotrans % ptr ( ics : ice , iblk )
!$acc data &
!$acc present(drag_wtr(:),q_acoef_wtr(:),q_bcoef_ice(:),q_air(:),evapo_wtr(:) &
!$acc ,qsat_lwtr(:),qsat_lice(:),drag_ice(:),evapo_ice(:),heat_tcoef(:) &
!$acc ,evapopot(:),evapotrans(:),q_bcoef_wtr(:),q_acoef_ice(:),fract_lice(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_6 = 1 , size ( q_air , 1 )
    q_air ( claw_induction_6 ) = q_acoef_wtr ( claw_induction_6 ) * qsat_lwtr&
     ( claw_induction_6 ) + q_bcoef_wtr ( claw_induction_6 )
    heat_tcoef ( claw_induction_6 ) = heat_transfer_coef ( drag_wtr (&
     claw_induction_6 ) , steplen )
    evapo_wtr ( claw_induction_6 ) = heat_tcoef ( claw_induction_6 ) * ( q_air&
     ( claw_induction_6 ) - qsat_lwtr ( claw_induction_6 ) )
    q_air ( claw_induction_6 ) = q_acoef_ice ( claw_induction_6 ) * qsat_lice&
     ( claw_induction_6 ) + q_bcoef_ice ( claw_induction_6 )
    heat_tcoef ( claw_induction_6 ) = heat_transfer_coef ( drag_ice (&
     claw_induction_6 ) , steplen )
    evapo_ice ( claw_induction_6 ) = heat_tcoef ( claw_induction_6 ) * ( q_air&
     ( claw_induction_6 ) - qsat_lice ( claw_induction_6 ) )
    evapopot ( claw_induction_6 ) = ( 1._wp - fract_lice ( claw_induction_6 )&
     ) * evapo_wtr ( claw_induction_6 ) + fract_lice ( claw_induction_6 ) *&
     evapo_ice ( claw_induction_6 )
    evapotrans ( claw_induction_6 ) = evapopot ( claw_induction_6 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
   IF ( tile % contains_soil .OR. tile % contains_glacier ) THEN
    q_acoef => a2l__mem % q_acoef % ptr ( ics : ice , iblk )
    q_bcoef => a2l__mem % q_bcoef % ptr ( ics : ice , iblk )
    drag_srf => a2l__mem % drag_srf % ptr ( ics : ice , iblk )
    qsat_star => seb__mem % qsat_star % ptr ( ics : ice , iblk )
    fact_q_air => turb__mem % fact_q_air % ptr ( ics : ice , iblk )
    fact_qsat_srf => turb__mem % fact_qsat_srf % ptr ( ics : ice , iblk )
    fact_qsat_trans_srf => turb__mem % fact_qsat_trans_srf % ptr ( ics : ice ,&
     iblk )
    evapotrans_lnd => hydro__mem % evapotrans_lnd % ptr ( ics : ice , iblk )
    IF ( tile % contains_vegetation ) THEN
     trans => hydro__mem % trans % ptr ( ics : ice , iblk )
    END IF
!$acc data &
!$acc present(q_bcoef(:),q_acoef(:),qsat_star(:),q_air(:),drag_srf(:) &
!$acc ,fact_q_air(:),qsat_srf_eff(:),evapotrans_lnd(:),heat_tcoef(:) &
!$acc ,evapopot(:),fact_qsat_srf(:),q_air_eff(:),evapotrans(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_7 = 1 , size ( q_air , 1 )
     q_air ( claw_induction_7 ) = q_acoef ( claw_induction_7 ) * qsat_star (&
      claw_induction_7 ) + q_bcoef ( claw_induction_7 )
     heat_tcoef ( claw_induction_7 ) = heat_transfer_coef ( drag_srf (&
      claw_induction_7 ) , steplen )
     q_air_eff ( claw_induction_7 ) = q_effective ( 0._wp , q_air (&
      claw_induction_7 ) , 1._wp , 0._wp )
     qsat_srf_eff ( claw_induction_7 ) = q_effective ( qsat_star (&
      claw_induction_7 ) , q_air ( claw_induction_7 ) , fact_qsat_srf (&
      claw_induction_7 ) , fact_q_air ( claw_induction_7 ) )
     evapotrans_lnd ( claw_induction_7 ) = heat_tcoef ( claw_induction_7 ) * (&
      q_air_eff ( claw_induction_7 ) - qsat_srf_eff ( claw_induction_7 ) )
     evapotrans ( claw_induction_7 ) = evapotrans_lnd ( claw_induction_7 )
     evapopot ( claw_induction_7 ) = heat_tcoef ( claw_induction_7 ) * ( q_air&
      ( claw_induction_7 ) - qsat_star ( claw_induction_7 ) )
    END DO
!$acc end parallel
!$acc end data
    IF ( tile % contains_vegetation ) THEN
!$acc data present(fact_qsat_trans_srf(:),trans(:),evapopot(:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_10 = 1 , size ( trans , 1 )
      trans ( claw_induction_10 ) = fact_qsat_trans_srf ( claw_induction_10 )&
       * evapopot ( claw_induction_10 )
     END DO
!$acc end parallel
!$acc end data
    END IF
   ELSE
    CALL finish ( trim ("mo_hydro_interface" //":update_evaporation" ) ,&
    "Called for invalid lct_type" )
   END IF
  END IF
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_evaporation" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_evaporation

 SUBROUTINE aggregate_evaporation ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_evaporation"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_evaporation" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapotrans , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_evaporation" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapopot , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_evaporation" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % trans , ics , ice ,&
   iblk ,"mo_hydro_interface" //":aggregate_evaporation" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapotrans_lnd ,&
   ics , ice , iblk ,"mo_hydro_interface" //":aggregate_evaporation" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapo_wtr , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_evaporation" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % evapo_ice , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_evaporation" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_evaporation" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_evaporation

 SUBROUTINE update_snow_and_ice_hydrology ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_snow_and_ice_hydrology"

  IF ( options % nc > 0 ) THEN
   CONTINUE
  END IF
  IF ( ( .NOT. tile % is_process_active ( hydro_ ) ) ) THEN
   RETURN
  END IF
  IF ( tile % contains_lake ) THEN
  ELSE
   IF ( tile % contains_land ) THEN
   ELSE
    CALL finish ( trim ("mo_hydro_interface" //&
    ":update_snow_and_ice_hydrology" ) ,"Called for invalid lct_type" )
   END IF
  END IF
 END SUBROUTINE update_snow_and_ice_hydrology

 SUBROUTINE aggregate_snow_and_ice_hydrology ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_snow_and_ice_hydrology"
  INTEGER :: iblk

  iblk = options % iblk
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":aggregate_snow_and_ice_hydrology" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":aggregate_snow_and_ice_hydrology" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_snow_and_ice_hydrology

 SUBROUTINE update_snow_and_skin_fraction ( tile , options )
  USE mo_hydro_process , ONLY: calc_wskin_fractions_lice ,&
   calc_wskin_fractions_veg , calc_wskin_fractions_bare

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  REAL ( KIND= wp ) :: w_skin_canopy_max ( 1 : options % nc )
  REAL ( KIND= wp ) :: w_skin_max ( 1 : options % nc )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_snow_and_skin_fraction"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: dtime
  REAL ( KIND= wp ) :: steplen
  REAL ( KIND= wp ) :: config_w_skin_max
INTEGER, SAVE :: task_counter = 1
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_hydro_config ) , POINTER :: hydro__conf
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow ( : )
  REAL ( KIND= wp ) , POINTER :: fract_water ( : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc_max ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_can ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_soil ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_lice ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow_can ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow_soil ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow_lice ( : )
  REAL ( KIND= wp ) , POINTER :: w_skin ( : )
  REAL ( KIND= wp ) , POINTER :: oro_stddev ( : )
  REAL ( KIND= wp ) , POINTER :: drag_srf ( : )
  REAL ( KIND= wp ) , POINTER :: t ( : )
  REAL ( KIND= wp ) , POINTER :: press_srf ( : )
  REAL ( KIND= wp ) , POINTER :: q_air ( : )
  INTEGER :: claw_induction_8
  INTEGER :: claw_induction_21
  INTEGER :: claw_induction_22
  INTEGER :: claw_induction_23
  INTEGER :: claw_induction_24
  INTEGER :: claw_induction_25
  INTEGER :: claw_induction_26
  INTEGER :: claw_induction_35

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  dtime = options % dtime
  steplen = options % steplen
  IF ( ( .NOT. tile % is_process_active ( hydro_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":update_snow_and_skin_fraction" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( hydro_ ) % p % config )
  TYPE IS ( t_hydro_config )
   hydro__conf => cc
  END SELECT
  config_w_skin_max = hydro__conf % w_skin_max
!$acc enter data pcreate(w_skin_max, w_skin_canopy_max)
!$acc data present(w_skin_max(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_35 = 1 , size ( w_skin_max , 1 )
   w_skin_max ( claw_induction_35 ) = config_w_skin_max
  END DO
!$acc end parallel
!$acc end data
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  IF ( tile % is_lake ) THEN
   w_snow_lice => hydro__mem % w_snow_lice % ptr ( ics : ice , iblk )
   fract_snow_lice => hydro__mem % fract_snow_lice % ptr ( ics : ice , iblk )
   CALL calc_wskin_fractions_lice ( w_snow_lice ( : ) , fract_snow_lice ( : )&
    )
!$acc exit data delete(w_skin_max, w_skin_canopy_max)
   RETURN
  END IF
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  IF ( tile % is_vegetation .OR. tile % is_land ) THEN
   SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
   TYPE IS ( t_pheno_memory )
    pheno__mem => mm
   END SELECT
  END IF
  drag_srf => a2l__mem % drag_srf % ptr ( ics : ice , iblk )
  q_air => a2l__mem % q_air % ptr ( ics : ice , iblk )
  press_srf => a2l__mem % press_srf % ptr ( ics : ice , iblk )
  t => seb__mem % t % ptr ( ics : ice , iblk )
  IF ( tile % is_bare .OR. tile % is_vegetation .OR. tile % is_land ) THEN
   oro_stddev => hydro__mem % oro_stddev % ptr ( ics : ice , iblk )
   fract_water => hydro__mem % fract_water % ptr ( ics : ice , iblk )
  END IF
  IF ( tile % is_vegetation .OR. tile % is_land ) THEN
   lai => pheno__mem % lai % ptr ( ics : ice , iblk )
   fract_fpc_max => pheno__mem % fract_fpc_max % ptr ( ics : ice , iblk )
   fract_snow_can => hydro__mem % fract_snow_can % ptr ( ics : ice , iblk )
  END IF
  IF ( tile % is_vegetation .OR. tile % is_land .OR. tile % is_bare ) THEN
   w_skin => hydro__mem % w_skin % ptr ( ics : ice , iblk )
  END IF
  fract_snow => hydro__mem % fract_snow % ptr ( ics : ice , iblk )
  w_snow_soil => hydro__mem % w_snow_soil % ptr ( ics : ice , iblk )
  fract_snow_soil => hydro__mem % fract_snow_soil % ptr ( ics : ice , iblk )
  IF ( tile % is_vegetation .OR. tile % is_land ) THEN
!$acc data present(lai(:),w_skin_canopy_max(:),fract_fpc_max(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_26 = 1 , size ( w_skin_canopy_max , 1 )
    w_skin_canopy_max ( claw_induction_26 ) = config_w_skin_max * lai (&
     claw_induction_26 ) * fract_fpc_max ( claw_induction_26 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(w_skin_canopy_max(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_25 = 1 , size ( w_skin_canopy_max , 1 )
    w_skin_canopy_max ( claw_induction_25 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( tile % contains_soil ) THEN
!$acc data present(w_skin_max(:),w_skin_canopy_max(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_24 = 1 , size ( w_skin_max , 1 )
    w_skin_max ( claw_induction_24 ) = config_w_skin_max + w_skin_canopy_max (&
     claw_induction_24 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(w_skin_max(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_23 = 1 , size ( w_skin_max , 1 )
    w_skin_max ( claw_induction_23 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( tile % is_glacier ) THEN
!$acc data present(fract_snow_soil(:),fract_snow(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_8 = 1 , size ( fract_snow , 1 )
    fract_snow ( claw_induction_8 ) = 1._wp
    fract_snow_soil ( claw_induction_8 ) = 1._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( tile % is_vegetation .OR. tile % is_land ) THEN
   w_snow_can => hydro__mem % w_snow_can % ptr ( ics : ice , iblk )
   CALL calc_wskin_fractions_veg ( dtime , steplen , w_skin_max ( : ) ,&
    oro_stddev ( : ) , t ( : ) , press_srf ( : ) , drag_srf ( : ) , q_air ( :&
    ) , w_skin ( : ) , w_snow_soil ( : ) , w_snow_can ( : ) , fract_snow_can (&
    : ) , fract_water ( : ) , fract_snow_soil ( : ) )
!$acc data present(fract_snow_soil(:),fract_snow(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_22 = 1 , size ( fract_snow , 1 )
    fract_snow ( claw_induction_22 ) = fract_snow_soil ( claw_induction_22 )
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( tile % is_bare ) THEN
   CALL calc_wskin_fractions_bare ( dtime , steplen , w_skin_max , oro_stddev&
    ( : ) , t ( : ) , press_srf ( : ) , drag_srf ( : ) , q_air ( : ) , w_skin&
    ( : ) , w_snow_soil ( : ) , fract_water ( : ) , fract_snow_soil ( : ) )
!$acc data present(fract_snow_soil(:),fract_snow(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_21 = 1 , size ( fract_snow , 1 )
    fract_snow ( claw_induction_21 ) = fract_snow_soil ( claw_induction_21 )
   END DO
!$acc end parallel
!$acc end data
  END IF
!$acc exit data delete(w_skin_max, w_skin_canopy_max)
 task_counter = task_counter + 1
 END SUBROUTINE update_snow_and_skin_fraction

 SUBROUTINE aggregate_snow_and_skin_fraction ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_snow_and_skin_fraction"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":aggregate_snow_and_skin_fraction" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % fract_water , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_snow_and_skin_fraction" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % fract_snow , ics ,&
   ice , iblk ,"mo_hydro_interface" //":aggregate_snow_and_skin_fraction" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % fract_snow_soil ,&
   ics , ice , iblk ,"mo_hydro_interface" //&
  ":aggregate_snow_and_skin_fraction" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % fract_snow_can ,&
   ics , ice , iblk ,"mo_hydro_interface" //&
  ":aggregate_snow_and_skin_fraction" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % fract_snow_lice ,&
   ics , ice , iblk ,"mo_hydro_interface" //&
  ":aggregate_snow_and_skin_fraction" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //&
   ":aggregate_snow_and_skin_fraction" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_snow_and_skin_fraction

 SUBROUTINE update_water_balance ( tile , options )
  USE mo_jsb_physical_constants , ONLY: rhoh2o

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_jsb_grid ) , POINTER :: grid
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  REAL ( KIND= wp ) , POINTER :: rain ( : )
  REAL ( KIND= wp ) , POINTER :: snow ( : )
  REAL ( KIND= wp ) , POINTER :: runoff ( : )
  REAL ( KIND= wp ) , POINTER :: drainage ( : )
  REAL ( KIND= wp ) , POINTER :: evapotrans ( : )
  REAL ( KIND= wp ) , POINTER :: w_skin ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow ( : )
  REAL ( KIND= wp ) , POINTER :: water_flux ( : )
  REAL ( KIND= wp ) , POINTER :: water_content ( : )
  REAL ( KIND= wp ) , POINTER :: water_budget ( : )
  REAL ( KIND= wp ) , POINTER :: w_soil_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_ice_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: area ( : )
  REAL ( KIND= wp ) :: tile_fract ( 1 : options % nc )
  REAL ( KIND= wp ) :: water_content_new ( 1 : options % nc )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: dtime
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":update_water_balance"
  INTEGER :: claw_induction_9
  INTEGER :: claw_induction_18
  INTEGER :: claw_induction_19
  INTEGER :: claw_induction_20
  INTEGER :: claw_induction_34

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  dtime = options % dtime
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_water_balance" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  grid => get_grid ( model % grid_id )
  area => grid % area ( ics : ice , iblk )
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  rain => a2l__mem % rain % ptr ( ics : ice , iblk )
  snow => a2l__mem % snow % ptr ( ics : ice , iblk )
  runoff => hydro__mem % runoff % ptr ( ics : ice , iblk )
  drainage => hydro__mem % drainage % ptr ( ics : ice , iblk )
  IF ( tile % contains_soil .AND. ( .NOT. tile % is_glacier ) ) THEN
   w_soil_sl => hydro__mem % w_soil_sl % ptr ( ics : ice , : , iblk )
   w_ice_sl => hydro__mem % w_ice_sl % ptr ( ics : ice , : , iblk )
   w_skin => hydro__mem % w_skin % ptr ( ics : ice , iblk )
   w_snow => hydro__mem % w_snow % ptr ( ics : ice , iblk )
  END IF
  evapotrans => hydro__mem % evapotrans % ptr ( ics : ice , iblk )
  water_flux => hydro__mem % water_flux % ptr ( ics : ice , iblk )
  water_content => hydro__mem % water_content % ptr ( ics : ice , iblk )
  water_budget => hydro__mem % water_budget % ptr ( ics : ice , iblk )
!$acc data create(tile_fract, water_content_new)
  CALL tile % get_fraction ( ics , ice , iblk , fract = tile_fract ( : ) )
!$acc update device(tile_fract)
!$acc data &
!$acc present(snow(:),tile_fract(:),drainage(:),rain(:),runoff(:),area(:) &
!$acc ,water_flux(:),evapotrans(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_9 = 1 , size ( water_flux , 1 )
   water_flux ( claw_induction_9 ) = rain ( claw_induction_9 ) + snow (&
    claw_induction_9 ) + evapotrans ( claw_induction_9 ) - runoff (&
    claw_induction_9 ) - drainage ( claw_induction_9 )
   water_flux ( claw_induction_9 ) = water_flux ( claw_induction_9 ) *&
    tile_fract ( claw_induction_9 ) * area ( claw_induction_9 ) / rhoh2o
  END DO
!$acc end parallel
!$acc end data
  IF ( tile % contains_soil .AND. ( .NOT. tile % is_glacier ) ) THEN
!$acc data &
!$acc present(tile_fract(:),w_ice_sl(:,:),area(:),w_soil_sl(:,:) &
!$acc ,water_content_new(:),w_skin(:),w_snow(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_20 = 1 , size ( water_content_new , 1 )
    water_content_new ( claw_induction_20 ) = ( w_skin ( claw_induction_20 ) +&
     w_snow ( claw_induction_20 ) + sum ( w_soil_sl ( claw_induction_20 , : )&
     ) + sum ( w_ice_sl ( claw_induction_20 , : ) ) ) * tile_fract (&
     claw_induction_20 ) * area ( claw_induction_20 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(water_content_new(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_19 = 1 , size ( water_content_new , 1 )
    water_content_new ( claw_induction_19 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( ( .NOT. is_time_experiment_start ( ) ) ) THEN
!$acc data &
!$acc present(water_budget(:),water_content(:),water_flux(:) &
!$acc ,water_content_new(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_18 = 1 , size ( water_budget , 1 )
    water_budget ( claw_induction_18 ) = water_content ( claw_induction_18 ) +&
     water_flux ( claw_induction_18 ) * dtime - water_content_new (&
     claw_induction_18 )
   END DO
!$acc end parallel
!$acc end data
  END IF
!$acc data present(water_content(:),water_content_new(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_34 = 1 , size ( water_content , 1 )
   water_content ( claw_induction_34 ) = water_content_new ( claw_induction_34&
    )
  END DO
!$acc end parallel
!$acc end data
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":update_water_balance" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_water_balance

 SUBROUTINE aggregate_water_balance ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_interface" //&
  ":aggregate_water_balance"
  INTEGER :: iblk

  iblk = options % iblk
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_water_balance" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  CALL update_water_balance ( tile , options )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_hydro_interface" //":aggregate_water_balance" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_water_balance

END MODULE mo_hydro_interface

