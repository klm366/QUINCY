MODULE mo_sse_interface
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message
 USE mo_jsb_control , ONLY: debug_on
 USE mo_jsb_time , ONLY: is_time_experiment_start
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_grid_class , ONLY: t_jsb_vgrid
 USE mo_jsb_grid , ONLY: get_vgrid
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract , t_jsb_aggregator
 USE mo_jsb_process_class , ONLY: t_jsb_process
 USE mo_jsb_task_class , ONLY: t_jsb_process_task , t_jsb_task_options
 USE mo_jsb_process_class , ONLY: sse_ , seb_ , hydro_
 USE mo_sse_config_class , ONLY: t_sse_config
 USE mo_hydro_config_class , ONLY: t_hydro_config
 USE mo_sse_memory_class , ONLY: t_sse_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_soil_and_snow_temp
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_soil_and_snow_temperature
  PROCEDURE , NOPASS :: aggregate => aggregate_soil_and_snow_temperature
 END TYPE tsk_soil_and_snow_temp
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_soil_and_snow_properties
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_soil_and_snow_properties
  PROCEDURE , NOPASS :: aggregate => aggregate_soil_and_snow_properties
 END TYPE tsk_soil_and_snow_properties
 PRIVATE :: create_task_soil_and_snow_temperature
 PRIVATE :: create_task_soil_and_snow_properties
 PRIVATE :: update_soil_and_snow_temperature
 PRIVATE :: aggregate_soil_and_snow_temperature
 PRIVATE :: update_soil_and_snow_properties
 PRIVATE :: aggregate_soil_and_snow_properties
 PUBLIC :: register_sse_tasks
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_sse_interface"
 INTERFACE tsk_soil_and_snow_temp
  PROCEDURE create_task_soil_and_snow_temperature
 END INTERFACE
 INTERFACE tsk_soil_and_snow_properties
  PROCEDURE create_task_soil_and_snow_properties
 END INTERFACE

CONTAINS
 FUNCTION create_task_soil_and_snow_temperature ( model_id )&
  RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_soil_and_snow_temp :: return_ptr )
  CALL return_ptr % construct ( name ="soil_and_snow_temperature" , process_id&
   = sse_ , owner_model_id = model_id )
 END FUNCTION create_task_soil_and_snow_temperature

 FUNCTION create_task_soil_and_snow_properties ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_soil_and_snow_properties :: return_ptr )
  CALL return_ptr % construct ( name ="soil_and_snow_properties" , process_id&
   = sse_ , owner_model_id = model_id )
 END FUNCTION create_task_soil_and_snow_properties

 SUBROUTINE register_sse_tasks ( process , model_id )

  CLASS ( t_jsb_process ) , INTENT(INOUT) :: process
  INTEGER , INTENT(IN) :: model_id

  CALL process % register_task ( tsk_soil_and_snow_temp ( model_id ) )
  CALL process % register_task ( tsk_soil_and_snow_properties ( model_id ) )
 END SUBROUTINE register_sse_tasks

 SUBROUTINE update_soil_and_snow_temperature ( tile , options )
  USE mo_sse_process , ONLY: calc_soil_temperature , calc_snow_temperature ,&
   calc_soiltemp_old
  USE mo_jsb_physical_constants , ONLY: alf , rhoh2o , dens_snow
  USE mo_sse_constants , ONLY: snow_depth_min

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: t_unfilt ( : )
  REAL ( KIND= wp ) , POINTER :: snowmelt ( : )
  REAL ( KIND= wp ) , POINTER :: t_soil ( : , : )
  REAL ( KIND= wp ) , POINTER :: t_soil_acoef ( : , : )
  REAL ( KIND= wp ) , POINTER :: t_soil_bcoef ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: heat_cond_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: t_snow ( : , : )
  REAL ( KIND= wp ) , POINTER :: t_snow_acoef ( : , : )
  REAL ( KIND= wp ) , POINTER :: t_snow_bcoef ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap_snow ( : )
  REAL ( KIND= wp ) , POINTER :: heat_cond_snow ( : )
  REAL ( KIND= wp ) , POINTER :: hcap_grnd ( : )
  REAL ( KIND= wp ) , POINTER :: hcap_grnd_old ( : )
  REAL ( KIND= wp ) , POINTER :: grnd_hflx ( : )
  REAL ( KIND= wp ) , POINTER :: grnd_hflx_old ( : )
  REAL ( KIND= wp ) , POINTER :: snow_depth_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: soil_depth_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_ice_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_freeze_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_ice_melt_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_porosity_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: matrix_pot_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: bclapp_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: thaw_depth ( : )
  REAL ( KIND= wp ) , POINTER :: heat_cond ( : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow_soil ( : )
  REAL ( KIND= wp ) , POINTER :: snow_soil_dens ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_soil ( : )
  REAL ( KIND= wp ) :: snow_depth ( 1 : options % nc )
  REAL ( KIND= wp ) :: snow_fract ( 1 : options % nc )
  REAL ( KIND= wp ) :: t_unfilt_corrected ( 1 : options % nc )
  REAL ( KIND= wp ) :: t_soil_top ( 1 : options % nc )
  REAL ( KIND= wp ) :: hcap_grnd_snow ( 1 : options % nc )
  REAL ( KIND= wp ) :: hcap_grnd_soil ( 1 : options % nc )
  REAL ( KIND= wp ) :: grnd_hflx_snow ( 1 : options % nc )
  REAL ( KIND= wp ) :: grnd_hflx_soil ( 1 : options % nc )
  LOGICAL :: ldglac ( 1 : options % nc )
  INTEGER :: itop ( 1 : options % nc )
  INTEGER :: itop_old ( 1 : options % nc )
  INTEGER :: isnow
  LOGICAL :: lstart
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: nc
  INTEGER :: i
  REAL ( KIND= wp ) :: dtime
  TYPE ( t_jsb_vgrid ) , POINTER :: soil_e
  TYPE ( t_jsb_vgrid ) , POINTER :: snow_e
  INTEGER :: nsoil
  INTEGER :: nsnow
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_interface" //&
  ":update_soil_and_snow_temperature"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_11
  INTEGER :: claw_induction_12
  INTEGER :: claw_induction_19
  INTEGER :: claw_induction_20
  INTEGER :: claw_induction_21
  INTEGER :: claw_induction_22
  INTEGER :: claw_induction_23
  INTEGER :: claw_induction_24
  INTEGER :: claw_induction_25
  INTEGER :: claw_induction_26

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  nc = options % nc
  dtime = options % dtime
  IF ( ( .NOT. tile % is_process_active ( sse_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_sse_interface" //&
   ":update_soil_and_snow_temperature" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( sse_ ) % p % config )
  TYPE IS ( t_sse_config )
   sse__conf => cc
  END SELECT
  soil_e => get_vgrid ("soil_depth_energy" )
  nsoil = soil_e % n_levels
  IF ( sse__conf % l_snow ) THEN
   snow_e => get_vgrid ("snow_depth_energy" )
   nsnow = snow_e % n_levels
  END IF
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  snowmelt => hydro__mem % snowmelt % ptr ( ics : ice , iblk )
  w_snow_soil => hydro__mem % w_snow_soil % ptr ( ics : ice , iblk )
  snow_soil_dens => hydro__mem % snow_soil_dens % ptr ( ics : ice , iblk )
  fract_snow_soil => hydro__mem % fract_snow_soil % ptr ( ics : ice , iblk )
  t_soil => sse__mem % t_soil % ptr ( ics : ice , : , iblk )
  t_soil_acoef => sse__mem % t_soil_acoef % ptr ( ics : ice , : , iblk )
  t_soil_bcoef => sse__mem % t_soil_bcoef % ptr ( ics : ice , : , iblk )
  vol_heat_cap_sl => sse__mem % vol_heat_cap_sl % ptr ( ics : ice , : , iblk )
  heat_cond_sl => sse__mem % heat_cond_sl % ptr ( ics : ice , : , iblk )
  vol_heat_cap_snow => sse__mem % vol_heat_cap_snow % ptr ( ics : ice , iblk )
  heat_cond_snow => sse__mem % heat_cond_snow % ptr ( ics : ice , iblk )
  hcap_grnd => sse__mem % hcap_grnd % ptr ( ics : ice , iblk )
  hcap_grnd_old => sse__mem % hcap_grnd_old % ptr ( ics : ice , iblk )
  grnd_hflx => sse__mem % grnd_hflx % ptr ( ics : ice , iblk )
  grnd_hflx_old => sse__mem % grnd_hflx_old % ptr ( ics : ice , iblk )
  IF ( sse__conf % l_snow ) THEN
   snow_depth_sl => sse__mem % snow_depth_sl % ptr ( ics : ice , : , iblk )
   t_snow => sse__mem % t_snow % ptr ( ics : ice , : , iblk )
   t_snow_acoef => sse__mem % t_snow_acoef % ptr ( ics : ice , : , iblk )
   t_snow_bcoef => sse__mem % t_snow_bcoef % ptr ( ics : ice , : , iblk )
  END IF
  IF ( ( .NOT. tile % is_glacier ) ) THEN
   soil_depth_sl => hydro__mem % soil_depth_sl % ptr ( ics : ice , : , iblk )
   vol_porosity_sl => hydro__mem % vol_porosity_sl % ptr ( ics : ice , : ,&
    iblk )
   matrix_pot_sl => hydro__mem % matrix_pot_sl % ptr ( ics : ice , : , iblk )
   bclapp_sl => hydro__mem % bclapp_sl % ptr ( ics : ice , : , iblk )
   w_soil_sl => hydro__mem % w_soil_sl % ptr ( ics : ice , : , iblk )
   w_ice_sl => hydro__mem % w_ice_sl % ptr ( ics : ice , : , iblk )
   w_soil_freeze_sl => hydro__mem % w_soil_freeze_sl % ptr ( ics : ice , : ,&
    iblk )
   w_ice_melt_sl => hydro__mem % w_ice_melt_sl % ptr ( ics : ice , : , iblk )
  END IF
  thaw_depth => sse__mem % thaw_depth % ptr ( ics : ice , iblk )
  heat_cond => sse__mem % heat_cond % ptr ( ics : ice , iblk )
  vol_heat_cap => sse__mem % vol_heat_cap % ptr ( ics : ice , iblk )
  t_unfilt => seb__mem % t_unfilt % ptr ( ics : ice , iblk )
  lstart = is_time_experiment_start ( )
!$acc enter data &
!$acc create(snow_depth,snow_fract,t_unfilt_corrected,t_soil_top &
!$acc ,hcap_grnd_snow,hcap_grnd_soil,grnd_hflx_snow,grnd_hflx_soil,ldglac &
!$acc ,itop,itop_old)
!$acc data present(grnd_hflx(:),hcap_grnd(:),grnd_hflx_old(:),hcap_grnd_old(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_0 = 1 , size ( grnd_hflx_old , 1 )
   grnd_hflx_old ( claw_induction_0 ) = grnd_hflx ( claw_induction_0 )
   hcap_grnd_old ( claw_induction_0 ) = hcap_grnd ( claw_induction_0 )
  END DO
!$acc end parallel
!$acc end data
  IF ( lstart ) THEN
!$acc data present(t_unfilt(:),t_unfilt_corrected(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_26 = 1 , size ( t_unfilt_corrected , 1 )
    t_unfilt_corrected ( claw_induction_26 ) = t_unfilt ( claw_induction_26 )
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(t_unfilt(:),hcap_grnd(:),t_unfilt_corrected(:),snowmelt(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_25 = 1 , size ( t_unfilt_corrected , 1 )
    t_unfilt_corrected ( claw_induction_25 ) = t_unfilt ( claw_induction_25 )&
     - snowmelt ( claw_induction_25 ) * dtime * alf / hcap_grnd (&
     claw_induction_25 )
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( ( .NOT. sse__conf % l_snow ) ) THEN
   ldglac ( : ) = spread ( tile % is_glacier , dim = 1 , ncopies = nc )
!$acc update device(ldglac)
   CALL calc_soiltemp_old ( nc , nsoil , soil_e % dz ( : ) , soil_e % levels (&
    : ) , dtime , lstart , t_unfilt_corrected ( : ) , w_snow_soil ( : ) ,&
    heat_cond ( : ) , vol_heat_cap ( : ) , t_soil_acoef ( : , : ) ,&
    t_soil_bcoef ( : , : ) , t_soil ( : , : ) , hcap_grnd ( : ) , grnd_hflx (&
    : ) , ldglac ( : ) )
!$acc exit data &
!$acc delete(snow_depth,snow_fract,t_unfilt_corrected,t_soil_top &
!$acc ,hcap_grnd_snow,hcap_grnd_soil,grnd_hflx_snow,grnd_hflx_soil,ldglac &
!$acc ,itop,itop_old)
   RETURN
  END IF
  IF ( ( .NOT. tile % is_glacier ) ) THEN
   IF ( sse__conf % l_dynsnow ) THEN
!$acc data present(w_snow_soil(:),snow_soil_dens(:),snow_depth(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_12 = 1 , size ( snow_depth , 1 )
     snow_depth ( claw_induction_12 ) = w_snow_soil ( claw_induction_12 ) *&
      rhoh2o / snow_soil_dens ( claw_induction_12 )
    END DO
!$acc end parallel
!$acc end data
   ELSE
!$acc data present(w_snow_soil(:),snow_depth(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_11 = 1 , size ( snow_depth , 1 )
     snow_depth ( claw_induction_11 ) = w_snow_soil ( claw_induction_11 ) *&
      rhoh2o / dens_snow
    END DO
!$acc end parallel
!$acc end data
   END IF
!$acc data present(itop,snow_depth_sl)
!$acc parallel loop gang vector
   DO i = 1 , nc , 1
    itop_old ( i ) = nsnow + 1 - count ( snow_depth_sl ( i , : ) > 0._wp )
   END DO
!$acc end data
!$acc data present(snow_depth_sl(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_23 = 1 , size ( snow_depth_sl , 1 )
    DO claw_induction_24 = 1 , size ( snow_depth_sl , 2 )
     snow_depth_sl ( claw_induction_23 , claw_induction_24 ) = 0._wp
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc data present(fract_snow_soil(:),snow_fract(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_22 = 1 , size ( snow_fract , 1 )
    snow_fract ( claw_induction_22 ) = fract_snow_soil ( claw_induction_22 )
   END DO
!$acc end parallel
!$acc end data
!$acc data present(snow_depth,snow_fract)
!$acc parallel loop gang vector
   DO i = 1 , size ( snow_depth , dim = 1 ) , 1
    IF ( snow_depth ( i ) > snow_depth_min .AND. snow_fract ( i ) > epsilon (&
     1._wp ) ) THEN
     snow_depth ( i ) = snow_depth ( i ) / snow_fract ( i )
    ELSE
     snow_depth ( i ) = 0._wp
     snow_fract ( i ) = 0._wp
    END IF
   END DO
!$acc end data
!$acc data present(snow_depth,snow_depth_sl)
!$acc parallel
!$acc loop seq
   DO isnow = nsnow , 1 , (-1)
!$acc loop vector
    DO i = 1 , size ( snow_depth , dim = 1 ) , 1
     IF ( snow_depth ( i ) > 0._wp ) THEN
      snow_depth_sl ( i , isnow ) = min ( snow_depth ( i ) , snow_e % dz (&
       isnow ) )
      snow_depth ( i ) = snow_depth ( i ) - snow_depth_sl ( i , isnow )
     END IF
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc data present(snow_depth_sl(:,nsnow),snow_depth(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_21 = 1 , size ( snow_depth_sl , 1 )
    snow_depth_sl ( claw_induction_21 , nsnow ) = snow_depth_sl (&
     claw_induction_21 , nsnow ) + snow_depth ( claw_induction_21 )
   END DO
!$acc end parallel
!$acc end data
   itop ( : ) = nsnow + 1 - count ( snow_depth_sl ( : , : ) > 0._wp , dim = 2&
    )
   CALL calc_snow_temperature ( nc , nsnow , soil_e % dz ( 1 : 2 ) , dtime ,&
    lstart , snow_depth_sl ( : , : ) , itop_old ( : ) , itop ( : ) ,&
    t_unfilt_corrected ( : ) , vol_heat_cap_snow ( : ) , heat_cond_snow ( : )&
    , vol_heat_cap_sl ( : , 1 ) , heat_cond_sl ( : , 1 ) , t_soil ( : , 1 ) ,&
    t_soil_acoef ( : , 1 ) , t_soil_bcoef ( : , 1 ) , t_snow ( : , : ) ,&
    t_snow_acoef ( : , : ) , t_snow_bcoef ( : , : ) , t_soil_top ( : ) ,&
    hcap_grnd_snow ( : ) , grnd_hflx_snow ( : ) )
   DO i = 1 , size ( snow_depth_sl , dim = 1 ) , 1
    IF ( snow_depth_sl ( i , nsnow ) > 0._wp ) THEN
     t_soil_top ( i ) = snow_fract ( i ) * t_soil_top ( i ) + ( 1._wp -&
      snow_fract ( i ) ) * t_unfilt_corrected ( i )
    ELSE
     t_soil_top ( i ) = t_unfilt_corrected ( i )
    END IF
   END DO
  ELSE
!$acc data present(snow_depth_sl(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_19 = 1 , size ( snow_depth_sl , 1 )
    DO claw_induction_20 = 1 , size ( snow_depth_sl , 2 )
     snow_depth_sl ( claw_induction_19 , claw_induction_20 ) = 0._wp
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc data &
!$acc present(w_snow_soil(:),t_soil_top(:),grnd_hflx_snow(:) &
!$acc ,t_unfilt_corrected(:),snow_fract(:),hcap_grnd_snow(:),snow_depth(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_1 = 1 , size ( snow_depth , 1 )
    snow_depth ( claw_induction_1 ) = w_snow_soil ( claw_induction_1 ) *&
     rhoh2o / dens_snow
    snow_fract ( claw_induction_1 ) = 0._wp
    t_soil_top ( claw_induction_1 ) = t_unfilt_corrected ( claw_induction_1 )
    grnd_hflx_snow ( claw_induction_1 ) = 0._wp
    hcap_grnd_snow ( claw_induction_1 ) = 0._wp
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( tile % is_glacier ) THEN
   CALL calc_soil_temperature ( nc , nsoil , soil_e % dz ( : ) , soil_e %&
    levels ( : ) , dtime , lstart , sse__conf % l_freeze , sse__conf %&
    l_supercool , t_soil_top ( : ) , vol_heat_cap_sl ( : , : ) , heat_cond_sl&
    ( : , : ) , t_soil ( : , : ) , t_soil_acoef ( : , : ) , t_soil_bcoef ( : ,&
    : ) , hcap_grnd_soil ( : ) , grnd_hflx_soil ( : ) , thaw_depth ( : ) )
  ELSE
   CALL calc_soil_temperature ( nc , nsoil , soil_e % dz ( : ) , soil_e %&
    levels ( : ) , dtime , lstart , sse__conf % l_freeze , sse__conf %&
    l_supercool , t_soil_top ( : ) , vol_heat_cap_sl ( : , : ) , heat_cond_sl&
    ( : , : ) , t_soil ( : , : ) , t_soil_acoef ( : , : ) , t_soil_bcoef ( : ,&
    : ) , hcap_grnd_soil ( : ) , grnd_hflx_soil ( : ) , thaw_depth ( : ) ,&
    soil_depth_sl ( : , : ) * vol_porosity_sl ( : , : ) , matrix_pot_sl ( : ,&
    : ) , bclapp_sl ( : , : ) , w_soil_sl ( : , : ) , w_ice_sl ( : , : ) ,&
    w_soil_freeze_sl ( : , : ) , w_ice_melt_sl ( : , : ) )
  END IF
!$acc data &
!$acc present(grnd_hflx,hcap_grnd,snow_fract,grnd_hflx_snow,grnd_hflx_soil &
!$acc ,hcap_grnd_snow,hcap_grnd_soil)
!$acc parallel loop gang vector
  DO i = 1 , size ( snow_depth_sl , dim = 1 ) , 1
   IF ( snow_depth_sl ( i , nsnow ) > 0._wp ) THEN
    grnd_hflx ( i ) = snow_fract ( i ) * grnd_hflx_snow ( i ) + ( 1._wp -&
     snow_fract ( i ) ) * grnd_hflx_soil ( i )
    hcap_grnd ( i ) = snow_fract ( i ) * hcap_grnd_snow ( i ) + ( 1._wp -&
     snow_fract ( i ) ) * hcap_grnd_soil ( i )
   ELSE
    grnd_hflx ( i ) = grnd_hflx_soil ( i )
    hcap_grnd ( i ) = hcap_grnd_soil ( i )
   END IF
  END DO
!$acc end data
!$acc exit data &
!$acc delete(snow_depth,snow_fract,t_unfilt_corrected,t_soil_top &
!$acc ,hcap_grnd_snow,hcap_grnd_soil,grnd_hflx_snow,grnd_hflx_soil,ldglac &
!$acc ,itop,itop_old)
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_sse_interface" //&
   ":update_soil_and_snow_temperature" ) ,"Finished." )
  END IF
 END SUBROUTINE update_soil_and_snow_temperature

 SUBROUTINE aggregate_soil_and_snow_temperature ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_interface" //&
  ":aggregate_soil_and_snow_temperature"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_sse_interface" //&
   ":aggregate_soil_and_snow_temperature" ) ,"Starting on tile " // trim (&
    tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( sse_ ) % p % config )
  TYPE IS ( t_sse_config )
   sse__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % t_soil_acoef , ics ,&
   ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % t_soil_bcoef , ics ,&
   ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % t_soil , ics , ice ,&
   iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % grnd_hflx , ics , ice&
   , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % grnd_hflx_old , ics ,&
   ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % hcap_grnd , ics , ice&
   , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % hcap_grnd_old , ics ,&
   ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  IF ( sse__conf % l_snow ) THEN
   CALL weighted_by_fract % aggregate ( tile , sse__mem % snow_depth_sl , ics&
    , ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature"&
    )
   CALL weighted_by_fract % aggregate ( tile , sse__mem % t_snow_acoef , ics ,&
    ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
   CALL weighted_by_fract % aggregate ( tile , sse__mem % t_snow_bcoef , ics ,&
    ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
   CALL weighted_by_fract % aggregate ( tile , sse__mem % t_snow , ics , ice ,&
    iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  END IF
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_soil_freeze_sl ,&
   ics , ice , iblk ,"mo_sse_interface" //&
  ":aggregate_soil_and_snow_temperature" )
  CALL weighted_by_fract % aggregate ( tile , hydro__mem % w_ice_melt_sl , ics&
   , ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % thaw_depth , ics ,&
   ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_temperature" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_sse_interface" //&
   ":aggregate_soil_and_snow_temperature" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_soil_and_snow_temperature

 SUBROUTINE update_soil_and_snow_properties ( tile , options )
  USE mo_jsb_physical_constants , ONLY: ci , dens_snow_min
  USE mo_sse_constants , ONLY: vol_hcap_ice , vol_hcap_snow , hcond_ice ,&
   hcond_snow , hcond_snow_min
  USE mo_sse_process , ONLY: calc_vol_heat_capacity ,&
   calc_thermal_conductivity

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_hydro_config ) , POINTER :: hydro__conf
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: snow_soil_dens ( : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap ( : )
  REAL ( KIND= wp ) , POINTER :: heat_cond ( : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap_snow ( : )
  REAL ( KIND= wp ) , POINTER :: heat_cond_snow ( : )
  REAL ( KIND= wp ) , POINTER :: soil_depth_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: fract_org_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: heat_cond_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: w_ice_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_porosity_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_field_cap_sl ( : , : )
  REAL ( KIND= wp ) , POINTER :: fract_water ( : , : )
  REAL ( KIND= wp ) , POINTER :: fract_ice ( : , : )
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_jsb_vgrid ) , POINTER :: soil_e
  REAL ( KIND= wp ) , ALLOCATABLE :: vol_heat_cap_tmp ( : , : )
  REAL ( KIND= wp ) , ALLOCATABLE :: soil_e_dz_tmp ( : , : )
  REAL ( KIND= wp ) , ALLOCATABLE :: heat_cond_tmp ( : , : )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: nc
  INTEGER :: i
  INTEGER :: j
  INTEGER :: nsoil
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_interface" //&
  ":update_soil_and_snow_properties"
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4
  INTEGER :: claw_induction_5
  INTEGER :: claw_induction_6
  INTEGER :: claw_induction_7
  INTEGER :: claw_induction_8
  INTEGER :: claw_induction_9
  INTEGER :: claw_induction_10
  INTEGER :: claw_induction_13
  INTEGER :: claw_induction_14
  INTEGER :: claw_induction_15
  INTEGER :: claw_induction_16
  INTEGER :: claw_induction_17
  INTEGER :: claw_induction_18

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  nc = options % nc
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_sse_interface" //&
   ":update_soil_and_snow_properties" ) ,"Starting on tile " // trim ( tile %&
    name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  soil_e => get_vgrid ("soil_depth_energy" )
  nsoil = soil_e % n_levels
  SELECT TYPE ( cc => model % processes ( sse_ ) % p % config )
  TYPE IS ( t_sse_config )
   sse__conf => cc
  END SELECT
  SELECT TYPE ( cc => model % processes ( hydro_ ) % p % config )
  TYPE IS ( t_hydro_config )
   hydro__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  vol_heat_cap_sl => sse__mem % vol_heat_cap_sl % ptr ( ics : ice , : , iblk )
  heat_cond_sl => sse__mem % heat_cond_sl % ptr ( ics : ice , : , iblk )
  vol_heat_cap_snow => sse__mem % vol_heat_cap_snow % ptr ( ics : ice , iblk )
  heat_cond_snow => sse__mem % heat_cond_snow % ptr ( ics : ice , iblk )
  IF ( tile % is_glacier ) THEN
   IF ( model % config % l_compat401 ) THEN
!$acc data present(vol_heat_cap_sl(:,:),heat_cond_sl(:,:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_2 = 1 , size ( vol_heat_cap_sl , 1 )
     DO claw_induction_3 = 1 , size ( vol_heat_cap_sl , 2 )
      vol_heat_cap_sl ( claw_induction_2 , claw_induction_3 ) = 2.09e+6_wp
      heat_cond_sl ( claw_induction_2 , claw_induction_3 ) = 2.508_wp
     END DO
    END DO
!$acc end parallel
!$acc end data
!$acc data present(vol_heat_cap_snow(:),heat_cond_snow(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_4 = 1 , size ( vol_heat_cap_snow , 1 )
     vol_heat_cap_snow ( claw_induction_4 ) = 634500._wp
     heat_cond_snow ( claw_induction_4 ) = 0.31_wp
    END DO
!$acc end parallel
!$acc end data
   ELSE
!$acc data present(vol_heat_cap_sl(:,:),heat_cond_sl(:,:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_5 = 1 , size ( vol_heat_cap_sl , 1 )
     DO claw_induction_6 = 1 , size ( vol_heat_cap_sl , 2 )
      vol_heat_cap_sl ( claw_induction_5 , claw_induction_6 ) = vol_hcap_ice
      heat_cond_sl ( claw_induction_5 , claw_induction_6 ) = hcond_ice
     END DO
    END DO
!$acc end parallel
!$acc end data
!$acc data present(vol_heat_cap_snow(:),heat_cond_snow(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_7 = 1 , size ( vol_heat_cap_snow , 1 )
     vol_heat_cap_snow ( claw_induction_7 ) = vol_hcap_snow
     heat_cond_snow ( claw_induction_7 ) = hcond_snow
    END DO
!$acc end parallel
!$acc end data
   END IF
   IF ( debug_on ( ) .AND. iblk == 1 ) THEN
    CALL message ( trim ("mo_sse_interface" //&
    ":update_soil_and_snow_properties" ) ,"Finished." )
   END IF
   RETURN
  END IF
  soil_depth_sl => hydro__mem % soil_depth_sl % ptr ( ics : ice , : , iblk )
  vol_heat_cap => sse__mem % vol_heat_cap % ptr ( ics : ice , iblk )
  heat_cond => sse__mem % heat_cond % ptr ( ics : ice , iblk )
  w_soil_sl => hydro__mem % w_soil_sl % ptr ( ics : ice , : , iblk )
  w_ice_sl => hydro__mem % w_ice_sl % ptr ( ics : ice , : , iblk )
  vol_porosity_sl => hydro__mem % vol_porosity_sl % ptr ( ics : ice , : , iblk&
   )
  vol_field_cap_sl => hydro__mem % vol_field_cap_sl % ptr ( ics : ice , : ,&
   iblk )
  IF ( hydro__conf % l_organic .AND. tile % contains_soil ) THEN
   fract_org_sl => hydro__mem % fract_org_sl % ptr ( ics : ice , : , iblk )
  ELSE
   ALLOCATE ( fract_org_sl ( nc , nsoil ) )
!$acc enter data pcreate(fract_org_sl)
!$acc data present(fract_org_sl(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_17 = 1 , size ( fract_org_sl , 1 )
    DO claw_induction_18 = 1 , size ( fract_org_sl , 2 )
     fract_org_sl ( claw_induction_17 , claw_induction_18 ) = 0._wp
    END DO
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( sse__conf % l_heat_cap_dyn ) THEN
   ALLOCATE ( fract_water ( nc , nsoil ) , fract_ice ( nc , nsoil ) )
!$acc enter data pcreate(fract_water, fract_ice)
!$acc data present(soil_depth_sl,fract_water,fract_ice,w_soil_sl,w_ice_sl)
!$acc parallel loop gang vector collapse(2)
   DO i = 1 , size ( soil_depth_sl , dim = 1 ) , 1
    DO j = 1 , size ( soil_depth_sl , dim = 2 ) , 1
     IF ( soil_depth_sl ( i , j ) > 0._wp ) THEN
      fract_water ( i , j ) = w_soil_sl ( i , j ) / soil_depth_sl ( i , j )
      fract_ice ( i , j ) = w_ice_sl ( i , j ) / soil_depth_sl ( i , j )
     ELSE
      fract_water ( i , j ) = 0._wp
      fract_ice ( i , j ) = 0._wp
     END IF
    END DO
   END DO
!$acc end data
   ALLOCATE ( vol_heat_cap_tmp ( nc , nsoil ) , soil_e_dz_tmp ( nc , nsoil ) )
!$acc enter data pcreate(vol_heat_cap_tmp,soil_e_dz_tmp)
!$acc update host(vol_heat_cap(:))
   vol_heat_cap_tmp ( : , : ) = spread ( vol_heat_cap ( : ) , dim = 2 ,&
    ncopies = nsoil )
   soil_e_dz_tmp ( : , : ) = spread ( soil_e % dz ( : ) , dim = 1 , ncopies =&
    nc )
!$acc update device(vol_heat_cap_tmp,soil_e_dz_tmp)
!$acc data &
!$acc present(vol_heat_cap_sl(:,:),vol_heat_cap_tmp(:,:),vol_porosity_sl(:,:) &
!$acc ,soil_depth_sl(:,:),fract_water(:,:),fract_ice(:,:),soil_e_dz_tmp(:,:) &
!$acc ,fract_org_sl(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_15 = 1 , size ( vol_heat_cap_sl , 1 )
    DO claw_induction_16 = 1 , size ( vol_heat_cap_sl , 2 )
     vol_heat_cap_sl ( claw_induction_15 , claw_induction_16 ) =&
      calc_vol_heat_capacity ( vol_heat_cap_tmp ( claw_induction_15 ,&
      claw_induction_16 ) , vol_porosity_sl ( claw_induction_15 ,&
      claw_induction_16 ) , fract_org_sl ( claw_induction_15 ,&
      claw_induction_16 ) , fract_water ( claw_induction_15 ,&
      claw_induction_16 ) , fract_ice ( claw_induction_15 , claw_induction_16&
      ) , soil_depth_sl ( claw_induction_15 , claw_induction_16 ) ,&
      soil_e_dz_tmp ( claw_induction_15 , claw_induction_16 ) )
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc exit data delete(fract_water,fract_ice,vol_heat_cap_tmp,soil_e_dz_tmp)
   DEALLOCATE ( fract_water , fract_ice , vol_heat_cap_tmp , soil_e_dz_tmp )
  ELSE
!$acc data present(vol_heat_cap_sl,vol_heat_cap)
!$acc parallel loop gang vector collapse(2)
   DO i = 1 , nc , 1
    DO j = 1 , nsoil , 1
     vol_heat_cap_sl ( i , j ) = vol_heat_cap ( i )
    END DO
   END DO
!$acc end data
  END IF
  IF ( sse__conf % l_heat_cond_dyn ) THEN
   ALLOCATE ( heat_cond_tmp ( nc , nsoil ) , soil_e_dz_tmp ( nc , nsoil ) )
!$acc enter data pcreate(heat_cond_tmp,soil_e_dz_tmp)
!$acc update host(heat_cond(:))
   heat_cond_tmp ( : , : ) = spread ( heat_cond ( : ) , dim = 2 , ncopies =&
    nsoil )
   soil_e_dz_tmp ( : , : ) = spread ( soil_e % dz ( : ) , dim = 1 , ncopies =&
    nc )
!$acc update device(heat_cond_tmp, soil_e_dz_tmp)
!$acc data &
!$acc present(heat_cond_tmp(:,:),vol_porosity_sl(:,:),soil_depth_sl(:,:) &
!$acc ,soil_e_dz_tmp(:,:),w_ice_sl(:,:),w_soil_sl(:,:),vol_field_cap_sl(:,:) &
!$acc ,fract_org_sl(:,:),heat_cond_sl(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_13 = 1 , size ( heat_cond_sl , 1 )
    DO claw_induction_14 = 1 , size ( heat_cond_sl , 2 )
     heat_cond_sl ( claw_induction_13 , claw_induction_14 ) =&
      calc_thermal_conductivity ( heat_cond_tmp ( claw_induction_13 ,&
      claw_induction_14 ) , w_soil_sl ( claw_induction_13 , claw_induction_14&
      ) , w_ice_sl ( claw_induction_13 , claw_induction_14 ) , vol_porosity_sl&
      ( claw_induction_13 , claw_induction_14 ) , vol_field_cap_sl (&
      claw_induction_13 , claw_induction_14 ) , fract_org_sl (&
      claw_induction_13 , claw_induction_14 ) , soil_depth_sl (&
      claw_induction_13 , claw_induction_14 ) , soil_e_dz_tmp (&
      claw_induction_13 , claw_induction_14 ) )
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc exit data delete(heat_cond_tmp,soil_e_dz_tmp)
   DEALLOCATE ( heat_cond_tmp , soil_e_dz_tmp )
  ELSE
!$acc data present(heat_cond_sl,heat_cond)
!$acc parallel loop gang vector collapse(2)
   DO i = 1 , nc , 1
    DO j = 1 , nsoil , 1
     heat_cond_sl ( i , j ) = heat_cond ( i )
    END DO
   END DO
!$acc end data
  END IF
  IF ( ( .NOT. ( hydro__conf % l_organic .AND. tile % contains_soil ) ) ) THEN
!$acc exit data delete(fract_org_sl)
   DEALLOCATE ( fract_org_sl )
  END IF
  IF ( model % config % l_compat401 ) THEN
!$acc data present(vol_heat_cap_snow(:),heat_cond_snow(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_8 = 1 , size ( vol_heat_cap_snow , 1 )
    vol_heat_cap_snow ( claw_induction_8 ) = 634500._wp
    heat_cond_snow ( claw_induction_8 ) = 0.31_wp
   END DO
!$acc end parallel
!$acc end data
  ELSE
   IF ( sse__conf % l_snow ) THEN
    IF ( sse__conf % l_dynsnow ) THEN
     snow_soil_dens => hydro__mem % snow_soil_dens % ptr ( ics : ice , iblk )
!$acc data present(snow_soil_dens,vol_heat_cap_snow,heat_cond_snow)
!$acc parallel loop gang vector
     DO i = 1 , size ( snow_soil_dens ) , 1
      IF ( snow_soil_dens ( i ) > dens_snow_min ) THEN
       vol_heat_cap_snow ( i ) = ci * snow_soil_dens ( i )
       heat_cond_snow ( i ) = 2.5e-6_wp * snow_soil_dens ( i ) ** ( 2._wp ) -&
        0.123e-3_wp * snow_soil_dens ( i ) + 0.024_wp
      ELSE
       heat_cond_snow ( i ) = hcond_snow_min
       vol_heat_cap_snow ( i ) = ci * dens_snow_min
      END IF
     END DO
!$acc end data
    ELSE
!$acc data present(vol_heat_cap_snow(:),heat_cond_snow(:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_9 = 1 , size ( vol_heat_cap_snow , 1 )
      vol_heat_cap_snow ( claw_induction_9 ) = vol_hcap_snow
      heat_cond_snow ( claw_induction_9 ) = hcond_snow
     END DO
!$acc end parallel
!$acc end data
    END IF
   ELSE
!$acc data present(vol_heat_cap_snow(:),heat_cond_snow(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_10 = 1 , size ( vol_heat_cap_snow , 1 )
     vol_heat_cap_snow ( claw_induction_10 ) = vol_hcap_snow
     heat_cond_snow ( claw_induction_10 ) = hcond_snow
    END DO
!$acc end parallel
!$acc end data
   END IF
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_sse_interface" //&
   ":update_soil_and_snow_properties" ) ,"Finished." )
  END IF
 END SUBROUTINE update_soil_and_snow_properties

 SUBROUTINE aggregate_soil_and_snow_properties ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_interface" //&
  ":aggregate_soil_and_snow_properties"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_sse_interface" //&
   ":aggregate_soil_and_snow_properties" ) ,"Starting on tile " // trim ( tile&
    % name ) //" ..." )
  END IF
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % vol_heat_cap_sl , ics&
   , ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_properties" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % heat_cond_sl , ics ,&
   ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_properties" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % vol_heat_cap_snow ,&
   ics , ice , iblk ,"mo_sse_interface" //&
  ":aggregate_soil_and_snow_properties" )
  CALL weighted_by_fract % aggregate ( tile , sse__mem % heat_cond_snow , ics&
   , ice , iblk ,"mo_sse_interface" //":aggregate_soil_and_snow_properties" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_sse_interface" //&
   ":aggregate_soil_and_snow_properties" ) ,"Finished." )
  END IF
 END SUBROUTINE aggregate_soil_and_snow_properties

END MODULE mo_sse_interface

