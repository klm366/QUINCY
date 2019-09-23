MODULE mo_seb_interface
 USE mo_jsb_control , ONLY: debug_on
 USE mo_jsb_time , ONLY: is_time_experiment_start
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message , finish
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_grid_class , ONLY: t_jsb_vgrid
 USE mo_jsb_grid , ONLY: get_vgrid
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract , t_jsb_aggregator
 USE mo_jsb_process_class , ONLY: t_jsb_process
 USE mo_jsb_task_class , ONLY: t_jsb_process_task , t_jsb_task_options
 USE mo_jsb_process_class , ONLY: seb_ , hydro_ , sse_
 USE mo_sse_config_class , ONLY: t_sse_config
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 USE mo_sse_memory_class , ONLY: t_sse_memory
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_surface_energy
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_surface_energy
  PROCEDURE , NOPASS :: aggregate => aggregate_surface_energy
 END TYPE tsk_surface_energy
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_snowmelt_correction
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_snowmelt_correction
  PROCEDURE , NOPASS :: aggregate => aggregate_snowmelt_correction
 END TYPE tsk_snowmelt_correction
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_asselin_filter
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_asselin
  PROCEDURE , NOPASS :: aggregate => aggregate_asselin
 END TYPE tsk_asselin_filter
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_surface_fluxes
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_surface_fluxes
  PROCEDURE , NOPASS :: aggregate => aggregate_surface_fluxes
 END TYPE tsk_surface_fluxes
 PRIVATE :: create_task_surface_energy
 PRIVATE :: create_task_snowmelt_correction
 PRIVATE :: create_task_asselin_filter
 PRIVATE :: create_task_surface_fluxes
 PRIVATE :: update_surface_energy
 PRIVATE :: aggregate_surface_energy
 PRIVATE :: update_snowmelt_correction
 PRIVATE :: aggregate_snowmelt_correction
 PRIVATE :: update_asselin
 PRIVATE :: aggregate_asselin
 PRIVATE :: update_surface_fluxes
 PRIVATE :: aggregate_surface_fluxes
 PUBLIC :: register_seb_tasks
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_seb_interface"
 INTERFACE tsk_surface_energy
  PROCEDURE create_task_surface_energy
 END INTERFACE
 INTERFACE tsk_snowmelt_correction
  PROCEDURE create_task_snowmelt_correction
 END INTERFACE
 INTERFACE tsk_asselin_filter
  PROCEDURE create_task_asselin_filter
 END INTERFACE
 INTERFACE tsk_surface_fluxes
  PROCEDURE create_task_surface_fluxes
 END INTERFACE

CONTAINS
 FUNCTION create_task_surface_energy ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_surface_energy :: return_ptr )
  CALL return_ptr % construct ( name ="surface_energy" , process_id = seb_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_surface_energy

 FUNCTION create_task_snowmelt_correction ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_snowmelt_correction :: return_ptr )
  CALL return_ptr % construct ( name ="snowmelt_correction" , process_id =&
   seb_ , owner_model_id = model_id )
 END FUNCTION create_task_snowmelt_correction

 FUNCTION create_task_asselin_filter ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_asselin_filter :: return_ptr )
  CALL return_ptr % construct ( name ="asselin_filter" , process_id = seb_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_asselin_filter

 FUNCTION create_task_surface_fluxes ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_surface_fluxes :: return_ptr )
  CALL return_ptr % construct ( name ="surface_fluxes" , process_id = seb_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_surface_fluxes

 SUBROUTINE register_seb_tasks ( this , model_id )

  CLASS ( t_jsb_process ) , INTENT(INOUT) :: this
  INTEGER , INTENT(IN) :: model_id

  CALL this % register_task ( tsk_surface_energy ( model_id ) )
  CALL this % register_task ( tsk_snowmelt_correction ( model_id ) )
  CALL this % register_task ( tsk_asselin_filter ( model_id ) )
  CALL this % register_task ( tsk_surface_fluxes ( model_id ) )
 END SUBROUTINE register_seb_tasks

 SUBROUTINE update_surface_energy ( tile , options )
  USE mo_seb_land , ONLY: update_surface_energy_land
  USE mo_seb_lake , ONLY: update_surface_energy_lake

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_interface" //&
  ":update_surface_energy"

  IF ( ( .NOT. tile % is_process_active ( seb_ ) ) ) THEN
   RETURN
  END IF
  IF ( tile % is_lake ) THEN
   CALL update_surface_energy_lake ( tile , options )
  ELSE
   IF ( tile % is_land ) THEN
    CALL update_surface_energy_land ( tile , options )
   ELSE
    CALL finish ( trim ("mo_seb_interface" //":update_surface_energy" ) ,&
    "Called for invalid lct_type " // trim ( tile % lcts ( 1 ) % name ) //&
    " on tile " // trim ( tile % name ) )
   END IF
  END IF
 END SUBROUTINE update_surface_energy

 SUBROUTINE aggregate_surface_energy ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_interface" //&
  ":aggregate_surface_energy"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":aggregate_surface_energy" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % t_unfilt , ics , ice&
   , iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % t_eff , ics , ice ,&
   iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % qsat_star , ics , ice&
   , iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % s_star , ics , ice ,&
   iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % forc_hflx , ics , ice&
   , iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % heat_cap , ics , ice&
   , iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % heat_cap_old , ics ,&
   ice , iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % t_lwtr , ics , ice ,&
   iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % t_lice , ics , ice ,&
   iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % s_wtr , ics , ice ,&
   iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % s_ice , ics , ice ,&
   iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % qsat_lwtr , ics , ice&
   , iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % qsat_lice , ics , ice&
   , iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % fract_lice , ics ,&
   ice , iblk ,"mo_seb_interface" //":aggregate_surface_energy" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":aggregate_surface_energy" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_surface_energy

 SUBROUTINE update_snowmelt_correction ( tile , options )
  USE mo_jsb_physical_constants , ONLY: alf

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: t_unfilt ( : )
  REAL ( KIND= wp ) , POINTER :: snowmelt ( : )
  REAL ( KIND= wp ) , POINTER :: hcap_grnd_old ( : )
  REAL ( KIND= wp ) , POINTER :: snow_depth_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_freeze_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_ice_melt_sl ( : , : )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: i
  REAL ( KIND= wp ) :: dtime
  REAL ( KIND= wp ) :: heat_cap ( 1 : options % nc )
  TYPE ( t_jsb_vgrid ) , POINTER :: soil_e
  TYPE ( t_jsb_vgrid ) , POINTER :: snow_e
  INTEGER :: nsnow
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_interface" //&
  ":update_snowmelt_correction"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1

  IF ( ( .NOT. tile % is_process_active ( seb_ ) ) .OR. ( .NOT. tile % is_land&
   ) ) THEN
   RETURN
  END IF
  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  dtime = options % dtime
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":update_snowmelt_correction" ) ,&
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
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  IF ( tile % contains_land ) THEN
   SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
   TYPE IS ( t_sse_memory )
    sse__mem => mm
   END SELECT
  END IF
  soil_e => get_vgrid ("soil_depth_energy" )
  IF ( sse__conf % l_snow ) THEN
   snow_e => get_vgrid ("snow_depth_energy" )
   nsnow = snow_e % n_levels
  END IF
  IF ( tile % contains_lake ) THEN
  ELSE
   IF ( tile % contains_land ) THEN
    t_unfilt => seb__mem % t_unfilt % ptr ( ics : ice , iblk )
    snowmelt => hydro__mem % snowmelt % ptr ( ics : ice , iblk )
    hcap_grnd_old => sse__mem % hcap_grnd_old % ptr ( ics : ice , iblk )
    IF ( sse__conf % l_snow ) THEN
     snow_depth_sl => sse__mem % snow_depth_sl % ptr ( ics : ice , : , iblk )
    END IF
    vol_heat_cap_sl => sse__mem % vol_heat_cap_sl % ptr ( ics : ice , : , iblk&
     )
    w_soil_freeze_sl => hydro__mem % w_soil_freeze_sl % ptr ( ics : ice , : ,&
     iblk )
    w_ice_melt_sl => hydro__mem % w_ice_melt_sl % ptr ( ics : ice , : , iblk )
    IF ( ( .NOT. is_time_experiment_start ( ) ) ) THEN
!$acc data create(heat_cap)
!$acc data present(t_unfilt(:),hcap_grnd_old(:),snowmelt(:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_0 = 1 , size ( t_unfilt , 1 )
      t_unfilt ( claw_induction_0 ) = t_unfilt ( claw_induction_0 ) - snowmelt&
       ( claw_induction_0 ) * dtime * alf / hcap_grnd_old ( claw_induction_0 )
     END DO
!$acc end parallel
!$acc end data
     IF ( sse__conf % l_freeze ) THEN
!$acc data present(heat_cap(:),vol_heat_cap_sl(:,1),soil_e%dz(1))
!$acc parallel
!$acc loop gang vector
      DO claw_induction_1 = 1 , size ( heat_cap , 1 )
       heat_cap ( claw_induction_1 ) = vol_heat_cap_sl ( claw_induction_1 , 1&
        ) * soil_e % dz ( 1 )
      END DO
!$acc end parallel
!$acc end data
      IF ( sse__conf % l_snow ) THEN
!$acc data present(t_unfilt,snow_depth_sl,heat_cap,w_ice_melt_sl,w_soil_freeze_sl)
!$acc parallel loop gang vector
       DO i = 1 , size ( t_unfilt , dim = 1 ) , 1
        IF ( snow_depth_sl ( i , nsnow ) < epsilon ( 1._wp ) .AND. heat_cap (&
         i ) > 0._wp ) THEN
        t_unfilt ( i ) = t_unfilt ( i ) - w_ice_melt_sl ( i , 1 ) * dtime *&
         alf / heat_cap ( i )
        t_unfilt ( i ) = t_unfilt ( i ) + w_soil_freeze_sl ( i , 1 ) * dtime *&
         alf / heat_cap ( i )
        END IF
       END DO
!$acc end data
      ELSE
!$acc data present(t_unfilt,heat_cap,w_ice_melt_sl,w_soil_freeze_sl)
!$acc parallel loop gang vector
       DO i = 1 , size ( t_unfilt , dim = 1 ) , 1
        IF ( heat_cap ( i ) > 0._wp ) THEN
        t_unfilt ( i ) = t_unfilt ( i ) - w_ice_melt_sl ( i , 1 ) * dtime *&
         alf / heat_cap ( i )
        t_unfilt ( i ) = t_unfilt ( i ) + w_soil_freeze_sl ( i , 1 ) * dtime *&
         alf / heat_cap ( i )
        END IF
       END DO
!$acc end data
      END IF
     END IF
!$acc end data
    END IF
   ELSE
    CALL finish ( trim ("mo_seb_interface" //":update_snowmelt_correction" ) ,&
    "Called for invalid lct_type" )
   END IF
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":update_snowmelt_correction" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_snowmelt_correction

 SUBROUTINE aggregate_snowmelt_correction ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_interface" //&
  ":aggregate_snowmelt_correction"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":aggregate_snowmelt_correction"&
    ) ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % t_unfilt , ics , ice&
   , iblk ,"mo_seb_interface" //":aggregate_snowmelt_correction" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":aggregate_snowmelt_correction"&
    ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_snowmelt_correction

 SUBROUTINE update_asselin ( tile , options )
  USE mo_seb_land , ONLY: update_asselin_land

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_interface" //&
  ":update_asselin"

  IF ( ( .NOT. tile % is_process_active ( seb_ ) ) ) THEN
   RETURN
  END IF
  IF ( tile % contains_lake ) THEN
  ELSE
   IF ( tile % contains_land ) THEN
    CALL update_asselin_land ( tile , options )
   ELSE
    CALL finish ( trim ("mo_seb_interface" //":update_asselin" ) ,&
    "Called for invalid lct_type" )
   END IF
  END IF
 END SUBROUTINE update_asselin

 SUBROUTINE aggregate_asselin ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_interface" //&
  ":aggregate_asselin"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":aggregate_asselin" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % t , ics , ice , iblk&
   ,"mo_seb_interface" //":aggregate_asselin" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % t_filt , ics , ice ,&
   iblk ,"mo_seb_interface" //":aggregate_asselin" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":aggregate_asselin" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_asselin

 SUBROUTINE update_surface_fluxes ( tile , options )
  USE mo_seb_land , ONLY: update_surface_fluxes_land
  USE mo_seb_lake , ONLY: update_surface_fluxes_lake

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_interface" //&
  ":update_surface_fluxes"

  IF ( ( .NOT. tile % is_process_active ( seb_ ) ) ) THEN
   RETURN
  END IF
  IF ( tile % contains_lake ) THEN
   CALL update_surface_fluxes_lake ( tile , options )
  ELSE
   IF ( tile % contains_land ) THEN
    CALL update_surface_fluxes_land ( tile , options )
   ELSE
    CALL finish ( trim ("mo_seb_interface" //":update_surface_fluxes" ) ,&
    "Called for invalid lct_type" )
   END IF
  END IF
 END SUBROUTINE update_surface_fluxes

 SUBROUTINE aggregate_surface_fluxes ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_interface" //&
  ":aggregate_surface_fluxes"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":aggregate_surface_fluxes" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % latent_hflx , ics ,&
   ice , iblk ,"mo_seb_interface" //":aggregate_surface_fluxes" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % sensible_hflx , ics ,&
   ice , iblk ,"mo_seb_interface" //":aggregate_surface_fluxes" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % latent_hflx_lnd , ics&
   , ice , iblk ,"mo_seb_interface" //":aggregate_surface_fluxes" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % sensible_hflx_lnd ,&
   ics , ice , iblk ,"mo_seb_interface" //":aggregate_surface_fluxes" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % latent_hflx_wtr , ics&
   , ice , iblk ,"mo_seb_interface" //":aggregate_surface_fluxes" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % sensible_hflx_wtr ,&
   ics , ice , iblk ,"mo_seb_interface" //":aggregate_surface_fluxes" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % latent_hflx_ice , ics&
   , ice , iblk ,"mo_seb_interface" //":aggregate_surface_fluxes" )
  CALL weighted_by_fract % aggregate ( tile , seb__mem % sensible_hflx_ice ,&
   ics , ice , iblk ,"mo_seb_interface" //":aggregate_surface_fluxes" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_interface" //":aggregate_surface_fluxes" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_surface_fluxes

END MODULE mo_seb_interface

