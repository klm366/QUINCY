MODULE mo_rad_interface
 USE mo_jsb_control , ONLY: debug_on
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message , finish
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract , t_jsb_aggregator
 USE mo_jsb_process_class , ONLY: t_jsb_process
 USE mo_jsb_task_class , ONLY: t_jsb_process_task , t_jsb_task_options
 USE mo_jsb_process_class , ONLY: seb_ , pheno_ , a2l_ , rad_ , hydro_ ,&
  assimi_ , carbon_
 USE mo_rad_config_class , ONLY: t_rad_config
 USE mo_assimi_config_class , ONLY: t_assimi_config
 USE mo_seb_config_class , ONLY: t_seb_config
 USE mo_a2l_memory_class , ONLY: t_a2l_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_pheno_memory_class , ONLY: t_pheno_memory
 USE mo_rad_memory_class , ONLY: t_rad_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 USE mo_carbon_memory_class , ONLY: t_carbon_memory
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_surface_radiation
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_surface_radiation
  PROCEDURE , NOPASS :: aggregate => aggregate_surface_radiation
 END TYPE tsk_surface_radiation
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_radiation_par
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_radiation_par
  PROCEDURE , NOPASS :: aggregate => aggregate_radiation_par
 END TYPE tsk_radiation_par
 TYPE , EXTENDS ( t_jsb_process_task ) :: tsk_albedo
 CONTAINS
  PROCEDURE , NOPASS :: integrate => update_albedo
  PROCEDURE , NOPASS :: aggregate => aggregate_albedo
 END TYPE tsk_albedo
 PRIVATE :: create_task_surface_radiation
 PRIVATE :: create_task_radiation_par
 PRIVATE :: create_task_albedo
 PRIVATE :: update_surface_radiation
 PRIVATE :: aggregate_surface_radiation
 PRIVATE :: update_radiation_par
 PRIVATE :: aggregate_radiation_par
 PRIVATE :: update_albedo
 PRIVATE :: aggregate_albedo
 PUBLIC :: register_rad_tasks
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_rad_interface"
 INTERFACE tsk_surface_radiation
  PROCEDURE create_task_surface_radiation
 END INTERFACE
 INTERFACE tsk_radiation_par
  PROCEDURE create_task_radiation_par
 END INTERFACE
 INTERFACE tsk_albedo
  PROCEDURE create_task_albedo
 END INTERFACE

CONTAINS
 FUNCTION create_task_surface_radiation ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_surface_radiation :: return_ptr )
  CALL return_ptr % construct ( name ="surface_radiation" , process_id = rad_&
   , owner_model_id = model_id )
 END FUNCTION create_task_surface_radiation

 FUNCTION create_task_radiation_par ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_radiation_par :: return_ptr )
  CALL return_ptr % construct ( name ="radiation_par" , process_id = rad_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_radiation_par

 FUNCTION create_task_albedo ( model_id ) RESULT(return_ptr)
  INTEGER , INTENT(IN) :: model_id
  CLASS ( t_jsb_process_task ) , POINTER :: return_ptr

  ALLOCATE ( tsk_albedo :: return_ptr )
  CALL return_ptr % construct ( name ="albedo" , process_id = rad_ ,&
   owner_model_id = model_id )
 END FUNCTION create_task_albedo

 SUBROUTINE register_rad_tasks ( this , model_id )

  CLASS ( t_jsb_process ) , INTENT(INOUT) :: this
  INTEGER , INTENT(IN) :: model_id

  CALL this % register_task ( tsk_surface_radiation ( model_id ) )
  CALL this % register_task ( tsk_radiation_par ( model_id ) )
  CALL this % register_task ( tsk_albedo ( model_id ) )
 END SUBROUTINE register_rad_tasks

 SUBROUTINE update_surface_radiation ( tile , options )
  USE mo_rad_process , ONLY: calc_radiation_surface_net

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_config ) , POINTER :: seb__conf
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  REAL ( KIND= wp ) , POINTER :: swvis_srf_net ( : )
  REAL ( KIND= wp ) , POINTER :: swnir_srf_net ( : )
  REAL ( KIND= wp ) , POINTER :: sw_srf_net ( : )
  REAL ( KIND= wp ) , POINTER :: lw_srf_net ( : )
  REAL ( KIND= wp ) , POINTER :: rad_srf_net ( : )
  REAL ( KIND= wp ) , POINTER :: lw_srf_down ( : )
  REAL ( KIND= wp ) , POINTER :: swvis_srf_down ( : )
  REAL ( KIND= wp ) , POINTER :: swnir_srf_down ( : )
  REAL ( KIND= wp ) , POINTER :: alb_vis_lnd ( : )
  REAL ( KIND= wp ) , POINTER :: alb_nir_lnd ( : )
  REAL ( KIND= wp ) , POINTER :: fract_lice ( : )
  REAL ( KIND= wp ) , POINTER :: rad_net_lwtr ( : )
  REAL ( KIND= wp ) , POINTER :: rad_net_lice ( : )
  REAL ( KIND= wp ) , POINTER :: albedo_lwtr ( : )
  REAL ( KIND= wp ) , POINTER :: albedo_lice ( : )
  REAL ( KIND= wp ) , POINTER :: t ( : )
  REAL ( KIND= wp ) , POINTER :: t_lwtr ( : )
  REAL ( KIND= wp ) , POINTER :: t_lice ( : )
  REAL ( KIND= wp ) :: lw_net_lwtr ( 1 : options % nc )
  REAL ( KIND= wp ) :: lw_net_lice ( 1 : options % nc )
  REAL ( KIND= wp ) :: sw_net_lwtr ( 1 : options % nc )
  REAL ( KIND= wp ) :: sw_net_lice ( 1 : options % nc )
  TYPE ( t_jsb_model ) , POINTER :: model
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_rad_interface" //&
  ":update_surface_radiation"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( ( .NOT. tile % is_process_active ( rad_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":update_surface_radiation" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( seb_ ) % p % config )
  TYPE IS ( t_seb_config )
   seb__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  swvis_srf_down => a2l__mem % swvis_srf_down % ptr ( ics : ice , iblk )
  swnir_srf_down => a2l__mem % swnir_srf_down % ptr ( ics : ice , iblk )
  lw_srf_down => a2l__mem % lw_srf_down % ptr ( ics : ice , iblk )
  sw_srf_net => rad__mem % sw_srf_net % ptr ( ics : ice , iblk )
  lw_srf_net => rad__mem % lw_srf_net % ptr ( ics : ice , iblk )
  rad_srf_net => rad__mem % rad_srf_net % ptr ( ics : ice , iblk )
  IF ( tile % is_lake ) THEN
   t_lwtr => seb__mem % t_lwtr % ptr ( ics : ice , iblk )
   albedo_lwtr => rad__mem % albedo_lwtr % ptr ( ics : ice , iblk )
   rad_net_lwtr => rad__mem % rad_net_lwtr % ptr ( ics : ice , iblk )
!$acc data pcreate(sw_net_lwtr, lw_net_lwtr, sw_net_lice, lw_net_lice)
   CALL calc_radiation_surface_net ( swvis_srf_down ( : ) , swnir_srf_down ( :&
    ) , albedo_lwtr ( : ) , albedo_lwtr ( : ) , lw_srf_down ( : ) , t_lwtr ( :&
    ) , rad_net_lwtr ( : ) , sw_net = sw_net_lwtr ( : ) , lw_net = lw_net_lwtr&
    ( : ) )
   IF ( seb__conf % l_ice_on_lakes ) THEN
    t_lice => seb__mem % t_lice % ptr ( ics : ice , iblk )
    fract_lice => seb__mem % fract_lice % ptr ( ics : ice , iblk )
    albedo_lice => rad__mem % albedo_lice % ptr ( ics : ice , iblk )
    rad_net_lice => rad__mem % rad_net_lice % ptr ( ics : ice , iblk )
    CALL calc_radiation_surface_net ( swvis_srf_down ( : ) , swnir_srf_down (&
     : ) , albedo_lice ( : ) , albedo_lice ( : ) , lw_srf_down ( : ) , t_lice&
     ( : ) , rad_net_lice ( : ) , sw_net = sw_net_lice ( : ) , lw_net =&
     lw_net_lice ( : ) )
!$acc data &
!$acc present(rad_net_lwtr(:),rad_net_lice(:),lw_net_lwtr(:),lw_net_lice(:) &
!$acc ,lw_srf_net(:),rad_srf_net(:),sw_srf_net(:),sw_net_lwtr(:) &
!$acc ,sw_net_lice(:),fract_lice(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_0 = 1 , size ( rad_srf_net , 1 )
     rad_srf_net ( claw_induction_0 ) = ( 1._wp - fract_lice (&
      claw_induction_0 ) ) * rad_net_lwtr ( claw_induction_0 ) + fract_lice (&
      claw_induction_0 ) * rad_net_lice ( claw_induction_0 )
     sw_srf_net ( claw_induction_0 ) = ( 1._wp - fract_lice ( claw_induction_0&
      ) ) * sw_net_lwtr ( claw_induction_0 ) + fract_lice ( claw_induction_0 )&
      * sw_net_lice ( claw_induction_0 )
     lw_srf_net ( claw_induction_0 ) = ( 1._wp - fract_lice ( claw_induction_0&
      ) ) * lw_net_lwtr ( claw_induction_0 ) + fract_lice ( claw_induction_0 )&
      * lw_net_lice ( claw_induction_0 )
    END DO
!$acc end parallel
!$acc end data
   ELSE
!$acc data &
!$acc present(rad_net_lwtr(:),lw_net_lwtr(:),lw_srf_net(:),rad_srf_net(:) &
!$acc ,sw_srf_net(:),sw_net_lwtr(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_1 = 1 , size ( rad_srf_net , 1 )
     rad_srf_net ( claw_induction_1 ) = rad_net_lwtr ( claw_induction_1 )
     sw_srf_net ( claw_induction_1 ) = sw_net_lwtr ( claw_induction_1 )
     lw_srf_net ( claw_induction_1 ) = lw_net_lwtr ( claw_induction_1 )
    END DO
!$acc end parallel
!$acc end data
   END IF
!$acc end data
  ELSE
   t => seb__mem % t % ptr ( ics : ice , iblk )
   alb_vis_lnd => rad__mem % alb_vis_lnd % ptr ( ics : ice , iblk )
   alb_nir_lnd => rad__mem % alb_nir_lnd % ptr ( ics : ice , iblk )
   swvis_srf_net => rad__mem % swvis_srf_net % ptr ( ics : ice , iblk )
   swnir_srf_net => rad__mem % swnir_srf_net % ptr ( ics : ice , iblk )
   CALL calc_radiation_surface_net ( swvis_srf_down ( : ) , swnir_srf_down ( :&
    ) , alb_vis_lnd ( : ) , alb_nir_lnd ( : ) , lw_srf_down ( : ) , t ( : ) ,&
    rad_srf_net ( : ) , swvis_srf_net ( : ) , swnir_srf_net ( : ) , sw_srf_net&
    ( : ) , lw_srf_net ( : ) )
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":update_surface_radiation" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_surface_radiation

 SUBROUTINE aggregate_surface_radiation ( tile , options )

  TYPE ( t_jsb_model ) , POINTER :: model
  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_config ) , POINTER :: seb__conf
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_rad_interface" //&
  ":aggregate_surface_radiation"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":aggregate_surface_radiation" )&
    ,"Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( seb_ ) % p % config )
  TYPE IS ( t_seb_config )
   seb__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % rad_srf_net , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_surface_radiation" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % lw_srf_net , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_surface_radiation" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % sw_srf_net , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_surface_radiation" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % swvis_srf_net , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_surface_radiation" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % swnir_srf_net , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_surface_radiation" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % rad_net_lwtr , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_surface_radiation" )
  IF ( seb__conf % l_ice_on_lakes ) THEN
   CALL weighted_by_fract % aggregate ( tile , rad__mem % rad_net_lice , ics ,&
    ice , iblk ,"mo_rad_interface" //":aggregate_surface_radiation" )
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":aggregate_surface_radiation" )&
    ,"Finished." )
  END IF
 END SUBROUTINE aggregate_surface_radiation

 SUBROUTINE update_radiation_par ( tile , options )
  USE mo_rad_process , ONLY: calc_par

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  INTEGER :: icanopy
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: dtime
  REAL ( KIND= wp ) :: canopy_bound_lai_delta
  REAL ( KIND= wp ) :: b4_layer_above ( 1 : options % nc )
  TYPE ( t_rad_config ) , POINTER :: rad__conf
  TYPE ( t_assimi_config ) , POINTER :: assimi__conf
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_jsb_model ) , POINTER :: model
  INTEGER :: ncanopy
  LOGICAL :: use_alb_veg_simple
  REAL ( KIND= wp ) , POINTER :: swpar_srf_down ( : )
  REAL ( KIND= wp ) , POINTER :: par ( : )
  REAL ( KIND= wp ) , POINTER :: par_down_mol ( : )
  REAL ( KIND= wp ) , POINTER :: par_down_mol_nacc ( : )
  REAL ( KIND= wp ) , POINTER :: soil_reflectivity_par ( : )
  REAL ( KIND= wp ) , POINTER :: alb_vis_soil ( : )
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: cos_zenith_angle ( : )
  REAL ( KIND= wp ) , POINTER :: fract_par_direct ( : )
  REAL ( KIND= wp ) , POINTER :: fract_par_diffuse ( : )
  REAL ( KIND= wp ) , POINTER :: canopy_bound_lai ( : )
  REAL ( KIND= wp ) , POINTER :: apar_nacc ( : )
  REAL ( KIND= wp ) , POINTER :: fapar ( : )
  REAL ( KIND= wp ) , POINTER :: fapar_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: apar_per_lai_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: lai_cl ( : , : )
  REAL ( KIND= wp ) , POINTER :: apar_per_lai_cl_nacc ( : , : )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_rad_interface" //&
  ":update_radiation_par"
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_7
  INTEGER :: claw_induction_8
  INTEGER :: claw_induction_10
  INTEGER :: claw_induction_11
  INTEGER :: claw_induction_12
  INTEGER :: claw_induction_13
  INTEGER :: claw_induction_14

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  dtime = options % dtime
  IF ( ( .NOT. tile % is_vegetation ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":update_radiation_par" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( rad_ ) % p % config )
  TYPE IS ( t_rad_config )
   rad__conf => cc
  END SELECT
  SELECT TYPE ( cc => model % processes ( assimi_ ) % p % config )
  TYPE IS ( t_assimi_config )
   assimi__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  IF ( ( .NOT. tile % is_process_active ( rad_ ) ) ) THEN
   RETURN
  END IF
  ncanopy = assimi__conf % ncanopy
  use_alb_veg_simple = rad__conf % use_alb_veg_simple
  canopy_bound_lai ( 0 : ncanopy ) => assimi__conf % canopy_bound_lai ( 0 :&
   ncanopy )
  cos_zenith_angle => a2l__mem % cos_zenith_angle % ptr ( ics : ice , iblk )
  fract_par_diffuse => a2l__mem % fract_par_diffuse % ptr ( ics : ice , iblk )
  swpar_srf_down => a2l__mem % swpar_srf_down % ptr ( ics : ice , iblk )
  lai => pheno__mem % lai % ptr ( ics : ice , iblk )
  fract_par_direct => rad__mem % fract_par_direct % ptr ( ics : ice , iblk )
  alb_vis_soil => rad__mem % alb_vis_soil % ptr ( ics : ice , iblk )
  par => rad__mem % par % ptr ( ics : ice , iblk )
  par_down_mol => rad__mem % par_down_mol % ptr ( ics : ice , iblk )
  par_down_mol_nacc => rad__mem % par_down_mol_nacc % ptr ( ics : ice , iblk )
  soil_reflectivity_par => rad__mem % soil_reflectivity_par % ptr ( ics : ice&
   , iblk )
  lai_cl => rad__mem % lai_cl % ptr ( ics : ice , : , iblk )
  fapar_cl => rad__mem % fapar_cl % ptr ( ics : ice , : , iblk )
  apar_per_lai_cl => rad__mem % apar_per_lai_cl % ptr ( ics : ice , : , iblk )
  apar_per_lai_cl_nacc => rad__mem % apar_per_lai_cl_nacc % ptr ( ics : ice ,&
   : , iblk )
  apar_nacc => rad__mem % apar_nacc % ptr ( ics : ice , iblk )
  fapar => rad__mem % fapar % ptr ( ics : ice , iblk )
!$acc data create(b4_layer_above)
!$acc data present(apar_per_lai_cl(:,:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_13 = 1 , size ( apar_per_lai_cl , 1 )
   DO claw_induction_14 = 1 , size ( apar_per_lai_cl , 2 )
    apar_per_lai_cl ( claw_induction_13 , claw_induction_14 ) = 0._wp
   END DO
  END DO
!$acc end parallel
!$acc end data
!$acc data &
!$acc present(soil_reflectivity_par(:),fract_par_diffuse(:) &
!$acc ,fract_par_direct(:),par_down_mol(:),b4_layer_above(:) &
!$acc ,swpar_srf_down(:),par(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_2 = 1 , size ( par_down_mol , 1 )
   par_down_mol ( claw_induction_2 ) = 0._wp
   soil_reflectivity_par ( claw_induction_2 ) = 0._wp
   par ( claw_induction_2 ) = swpar_srf_down ( claw_induction_2 )
   fract_par_direct ( claw_induction_2 ) = 1._wp - fract_par_diffuse (&
    claw_induction_2 )
   b4_layer_above ( claw_induction_2 ) = 7777777._wp
  END DO
!$acc end parallel
!$acc end data
  DO icanopy = 1 , ncanopy , 1
   canopy_bound_lai_delta = canopy_bound_lai ( icanopy ) - canopy_bound_lai (&
    icanopy - 1 )
   CALL calc_par ( ncanopy , icanopy , use_alb_veg_simple , par ( : ) ,&
    alb_vis_soil ( : ) , lai ( : ) , cos_zenith_angle ( : ) , fract_par_direct&
    ( : ) , canopy_bound_lai ( icanopy ) , canopy_bound_lai_delta ,&
    b4_layer_above ( : ) , par_down_mol ( : ) , soil_reflectivity_par ( : ) ,&
    fapar_cl ( : , icanopy ) , lai_cl ( : , icanopy ) , apar_per_lai_cl ( : ,&
    icanopy ) )
  END DO
!$acc data present(par_down_mol(:),par_down_mol_nacc(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_12 = 1 , size ( par_down_mol_nacc , 1 )
   par_down_mol_nacc ( claw_induction_12 ) = par_down_mol_nacc (&
    claw_induction_12 ) + par_down_mol ( claw_induction_12 ) * dtime
  END DO
!$acc end parallel
!$acc end data
!$acc data present(apar_per_lai_cl(:,:),apar_per_lai_cl_nacc(:,:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_10 = 1 , size ( apar_per_lai_cl_nacc , 1 )
   DO claw_induction_11 = 1 , size ( apar_per_lai_cl_nacc , 2 )
    apar_per_lai_cl_nacc ( claw_induction_10 , claw_induction_11 ) =&
     apar_per_lai_cl_nacc ( claw_induction_10 , claw_induction_11 ) +&
     apar_per_lai_cl ( claw_induction_10 , claw_induction_11 ) * dtime
   END DO
  END DO
!$acc end parallel
!$acc end data
  DO icanopy = 1 , ncanopy , 1
!$acc data present(par_down_mol(:),apar_nacc(:),fapar_cl(:,icanopy))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_8 = 1 , size ( apar_nacc , 1 )
    apar_nacc ( claw_induction_8 ) = apar_nacc ( claw_induction_8 ) +&
     par_down_mol ( claw_induction_8 ) * fapar_cl ( claw_induction_8 , icanopy&
     ) * dtime
   END DO
!$acc end parallel
!$acc end data
  END DO
  DO icanopy = 1 , ncanopy , 1
!$acc data present(fapar(:),fapar_cl(:,icanopy))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_7 = 1 , size ( fapar , 1 )
    fapar ( claw_induction_7 ) = fapar ( claw_induction_7 ) + fapar_cl (&
     claw_induction_7 , icanopy )
   END DO
!$acc end parallel
!$acc end data
  END DO
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":update_radiation_par" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_radiation_par

 SUBROUTINE aggregate_radiation_par ( tile , options )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_rad_interface" //&
  ":aggregate_radiation_par"
  INTEGER :: iblk

  iblk = options % iblk
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":aggregate_radiation_par" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":aggregate_radiation_par" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_radiation_par

 SUBROUTINE update_albedo ( tile , options )
  USE mo_rad_process , ONLY: get_surface_albedo_simple , calc_soil_albedo ,&
   calc_alb_lwtr , calc_alb_lice , calc_glacier_albedo ,&
   calc_sky_view_fractions , calc_snow_albedo , merge_albedos_of_vegtile ,&
   has_minimal_radiation
  USE mo_jsb_time , ONLY: is_time_experiment_start , is_time_ltrig_rad_m1
  USE mo_jsb_control , ONLY: jsbach_runs_standalone

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_jsb_model ) , POINTER :: model
  INTEGER :: iblk
  INTEGER :: nc
  INTEGER :: ics
  INTEGER :: ice
  INTEGER :: i
  LOGICAL , ALLOCATABLE :: l_day ( : )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_rad_interface" //&
  ":update_albedo"
  TYPE ( t_seb_config ) , POINTER :: seb__conf
  TYPE ( t_rad_config ) , POINTER :: rad__conf
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  TYPE ( t_carbon_memory ) , POINTER :: carbon__mem
  CHARACTER ( LEN= : ) , ALLOCATABLE :: use_alb_soil_organic_c
  INTEGER :: use_alb_soil_organic_int
  REAL ( KIND= wp ) :: albedo_age_weight
  REAL ( KIND= wp ) :: albedocanopysnow
  LOGICAL :: use_alb_veg_simple
  LOGICAL :: use_alb_soil_scheme
  LOGICAL :: use_alb_soil_litter
  LOGICAL :: use_alb_mineralsoil_const
  REAL ( KIND= wp ) , POINTER :: cos_zenith_angle ( : )
  REAL ( KIND= wp ) , POINTER :: t_unfilt ( : )
  REAL ( KIND= wp ) , POINTER :: alb_vis ( : )
  REAL ( KIND= wp ) , POINTER :: alb_nir ( : )
  REAL ( KIND= wp ) , POINTER :: alb_vis_lnd ( : )
  REAL ( KIND= wp ) , POINTER :: alb_nir_lnd ( : )
  REAL ( KIND= wp ) , POINTER :: alb_vis_can ( : )
  REAL ( KIND= wp ) , POINTER :: alb_nir_can ( : )
  REAL ( KIND= wp ) , POINTER :: alb_vis_soil ( : )
  REAL ( KIND= wp ) , POINTER :: alb_nir_soil ( : )
  REAL ( KIND= wp ) , POINTER :: alb_vis_snow ( : )
  REAL ( KIND= wp ) , POINTER :: alb_nir_snow ( : )
  REAL ( KIND= wp ) , POINTER :: alb_background ( : )
  REAL ( KIND= wp ) , POINTER :: alb_vis_mineralsoil ( : )
  REAL ( KIND= wp ) , POINTER :: alb_nir_mineralsoil ( : )
  REAL ( KIND= wp ) , POINTER :: sky_view_fract ( : )
  REAL ( KIND= wp ) , POINTER :: sky_view_fract_stem ( : )
  REAL ( KIND= wp ) , POINTER :: snow_age ( : )
  REAL ( KIND= wp ) , POINTER :: lai ( : )
  REAL ( KIND= wp ) , POINTER :: fract_fpc_max ( : )
  REAL ( KIND= wp ) , POINTER :: fract_forest ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_can ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_soil ( : )
  REAL ( KIND= wp ) , POINTER :: c_ag_sum_1 ( : )
  REAL ( KIND= wp ) , POINTER :: c_bg_sum ( : )
  REAL ( KIND= wp ) , POINTER :: swvis_srf_down ( : )
  REAL ( KIND= wp ) , POINTER :: swnir_srf_down ( : )
  REAL ( KIND= wp ) , POINTER :: albedo_lwtr ( : )
  REAL ( KIND= wp ) , POINTER :: fract_lice ( : )
  REAL ( KIND= wp ) , POINTER :: t_lice ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow_lice ( : )
  REAL ( KIND= wp ) , POINTER :: albedo_lice ( : )
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4
  INTEGER :: claw_induction_5
  INTEGER :: claw_induction_6
  INTEGER :: claw_induction_9

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  nc = options % nc
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":update_albedo" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  IF ( ( .NOT. is_time_experiment_start ( ) ) .AND. ( .NOT.&
   is_time_ltrig_rad_m1 ( model % id ) ) .AND. ( .NOT. jsbach_runs_standalone&
   ( ) ) ) THEN
   RETURN
  END IF
  SELECT TYPE ( cc => model % processes ( seb_ ) % p % config )
  TYPE IS ( t_seb_config )
   seb__conf => cc
  END SELECT
  SELECT TYPE ( cc => model % processes ( rad_ ) % p % config )
  TYPE IS ( t_rad_config )
   rad__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
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
  END IF
  use_alb_veg_simple = rad__conf % use_alb_veg_simple
  use_alb_soil_scheme = rad__conf % use_alb_soil_scheme
  use_alb_mineralsoil_const = rad__conf % use_alb_mineralsoil_const
  use_alb_soil_litter = rad__conf % use_alb_soil_litter
  use_alb_soil_organic_c = rad__conf % use_alb_soil_organic_c
  albedo_age_weight = rad__conf % albedo_age_weight
  albedocanopysnow = rad__conf % albedocanopysnow
  IF ( use_alb_soil_organic_c .eq."linear" ) THEN
   use_alb_soil_organic_int = 1
  ELSE
   IF ( use_alb_soil_organic_c .eq."log" ) THEN
    use_alb_soil_organic_int = 2
   END IF
  END IF
  c_bg_sum => null ( )
  IF ( tile % is_bare .OR. tile % is_vegetation ) THEN
   IF ( trim ( use_alb_soil_organic_c ) /="" .AND. ( .NOT. tile %&
    is_process_active ( carbon_ ) ) ) THEN
    CALL finish ( trim ("mo_rad_interface" //":update_albedo" ) ,&
    "Carbon process must be active for use_alb_soil_organic_C" )
   END IF
   IF ( trim ( use_alb_soil_organic_c ) /="" ) THEN
    SELECT TYPE ( mm => tile % mem ( carbon_ ) % p )
    TYPE IS ( t_carbon_memory )
     carbon__mem => mm
    END SELECT
    c_bg_sum => carbon__mem % c_bg_sum % ptr ( ics : ice , iblk )
   ELSE
    ALLOCATE ( c_bg_sum ( nc ) )
!$acc enter data create(c_bg_sum)
!$acc data present(c_bg_sum(:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_6 = 1 , size ( c_bg_sum , 1 )
     c_bg_sum ( claw_induction_6 ) = 0._wp
    END DO
!$acc end parallel
!$acc end data
   END IF
  END IF
  cos_zenith_angle => a2l__mem % cos_zenith_angle % ptr ( ics : ice , iblk )
  swvis_srf_down => a2l__mem % swvis_srf_down % ptr ( ics : ice , iblk )
  swnir_srf_down => a2l__mem % swnir_srf_down % ptr ( ics : ice , iblk )
  t_unfilt => seb__mem % t_unfilt % ptr ( ics : ice , iblk )
  alb_background => rad__mem % alb_background % ptr ( ics : ice , iblk )
  alb_vis => rad__mem % alb_vis % ptr ( ics : ice , iblk )
  alb_nir => rad__mem % alb_nir % ptr ( ics : ice , iblk )
  alb_vis_snow => rad__mem % alb_vis_snow % ptr ( ics : ice , iblk )
  alb_nir_snow => rad__mem % alb_nir_snow % ptr ( ics : ice , iblk )
  IF ( ( .NOT. tile % is_lake ) ) THEN
   alb_vis_lnd => rad__mem % alb_vis_lnd % ptr ( ics : ice , iblk )
   alb_nir_lnd => rad__mem % alb_nir_lnd % ptr ( ics : ice , iblk )
  END IF
  ALLOCATE ( l_day ( nc ) )
!$acc enter data create(l_day)
!$acc data present(swvis_srf_down(:),l_day(:),swnir_srf_down(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_9 = 1 , size ( l_day , 1 )
   l_day ( claw_induction_9 ) = has_minimal_radiation ( swvis_srf_down (&
    claw_induction_9 ) , swnir_srf_down ( claw_induction_9 ) )
  END DO
!$acc end parallel
!$acc end data
  IF ( tile % is_bare .OR. tile % is_vegetation ) THEN
   fract_snow_soil => hydro__mem % fract_snow_soil % ptr ( ics : ice , iblk )
   snow_age => rad__mem % snow_age % ptr ( ics : ice , iblk )
   CALL calc_snow_albedo ( l_day ( : ) , fract_snow_soil ( : ) , snow_age ( :&
    ) , cos_zenith_angle ( : ) , albedo_age_weight , t_unfilt ( : ) ,&
    alb_vis_snow ( : ) , alb_nir_snow ( : ) )
  END IF
  IF ( tile % is_lake ) THEN
   t_lice => seb__mem % t_lice % ptr ( ics : ice , iblk )
   albedo_lwtr => rad__mem % albedo_lwtr % ptr ( ics : ice , iblk )
   CALL calc_alb_lwtr ( l_day ( : ) , albedo_lwtr ( : ) )
   IF ( seb__conf % l_ice_on_lakes ) THEN
    fract_lice => seb__mem % fract_lice % ptr ( ics : ice , iblk )
    w_snow_lice => hydro__mem % w_snow_lice % ptr ( ics : ice , iblk )
    albedo_lice => rad__mem % albedo_lice % ptr ( ics : ice , iblk )
    CALL calc_alb_lice ( l_day ( : ) , t_lice ( : ) , w_snow_lice ( : ) ,&
     albedo_lice ( : ) )
!$acc data present(l_day,alb_vis,alb_nir,fract_lice,albedo_lwtr,albedo_lice)
!$acc parallel loop gang vector
    DO i = 1 , size ( l_day ) , 1
     IF ( l_day ( i ) ) THEN
      alb_vis ( i ) = ( 1._wp - fract_lice ( i ) ) * albedo_lwtr ( i ) +&
       fract_lice ( i ) * albedo_lice ( i )
      alb_nir ( i ) = alb_vis ( i )
     END IF
    END DO
!$acc end data
   ELSE
!$acc data present(l_day,alb_vis,alb_nir,fract_lice,albedo_lwtr,albedo_lice)
!$acc parallel loop gang vector
    DO i = 1 , size ( l_day ) , 1
     IF ( l_day ( i ) ) THEN
      alb_vis ( i ) = albedo_lwtr ( i )
      alb_nir ( i ) = alb_vis ( i )
     END IF
    END DO
!$acc end data
   END IF
  ELSE
   IF ( tile % is_glacier ) THEN
    CALL calc_glacier_albedo ( l_day ( : ) , t_unfilt ( : ) , alb_vis_lnd ( :&
     ) , alb_nir_lnd ( : ) )
   ELSE
    IF ( tile % is_bare ) THEN
     alb_vis_mineralsoil => rad__mem % alb_vis_mineralsoil % ptr ( ics : ice ,&
      iblk )
     alb_nir_mineralsoil => rad__mem % alb_nir_mineralsoil % ptr ( ics : ice ,&
      iblk )
     alb_vis_soil => rad__mem % alb_vis_soil % ptr ( ics : ice , iblk )
     alb_nir_soil => rad__mem % alb_nir_soil % ptr ( ics : ice , iblk )
     IF ( use_alb_soil_scheme ) THEN
      IF ( tile % lcts ( 1 ) % lib_id == 0 ) THEN
       CALL finish ( trim ("mo_rad_interface" //":update_albedo" ) ,&
       "lctlib not available on tile" )
      END IF
      c_ag_sum_1 => carbon__mem % c_ag_sum_1 % ptr ( ics : ice , iblk )
      CALL calc_soil_albedo ( l_day ( : ) , use_alb_mineralsoil_const ,&
       use_alb_soil_organic_int , use_alb_soil_litter , c_ag_sum_1 ( : ) ,&
       c_bg_sum ( : ) , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
       specificleafarea_c , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
       albedolittervis , alb_background ( : ) , alb_vis_mineralsoil ( : ) ,&
       alb_nir_mineralsoil ( : ) , alb_vis_soil ( : ) , alb_nir_soil ( : ) )
     ELSE
      IF ( debug_on ( ) .AND. iblk == 1 ) THEN
       CALL message ( trim ("mo_rad_interface" //":update_albedo" ) ,&
       "Soil albedo scheme not used for bare soil tiles but taken from bc_land&
       &_phy.nc file" )
      END IF
     END IF
!$acc data &
!$acc present(alb_vis_soil(:),alb_vis_lnd(:),alb_vis_snow(:),alb_nir_lnd(:) &
!$acc ,fract_snow_soil(:),alb_nir_soil(:),alb_nir_snow(:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_3 = 1 , size ( alb_vis_lnd , 1 )
      alb_vis_lnd ( claw_induction_3 ) = fract_snow_soil ( claw_induction_3 )&
       * alb_vis_snow ( claw_induction_3 ) + ( 1.0_wp - fract_snow_soil (&
       claw_induction_3 ) ) * alb_vis_soil ( claw_induction_3 )
      alb_nir_lnd ( claw_induction_3 ) = fract_snow_soil ( claw_induction_3 )&
       * alb_nir_snow ( claw_induction_3 ) + ( 1.0_wp - fract_snow_soil (&
       claw_induction_3 ) ) * alb_nir_soil ( claw_induction_3 )
     END DO
!$acc end parallel
!$acc end data
    ELSE
     IF ( tile % is_vegetation ) THEN
      lai => pheno__mem % lai % ptr ( ics : ice , iblk )
      fract_fpc_max => pheno__mem % fract_fpc_max % ptr ( ics : ice , iblk )
      fract_forest => pheno__mem % fract_forest % ptr ( ics : ice , iblk )
      fract_snow_can => hydro__mem % fract_snow_can % ptr ( ics : ice , iblk )
      fract_snow_soil => hydro__mem % fract_snow_soil % ptr ( ics : ice , iblk&
       )
      alb_vis_mineralsoil => rad__mem % alb_vis_mineralsoil % ptr ( ics : ice&
       , iblk )
      alb_nir_mineralsoil => rad__mem % alb_nir_mineralsoil % ptr ( ics : ice&
       , iblk )
      alb_vis_can => rad__mem % alb_vis_can % ptr ( ics : ice , iblk )
      alb_nir_can => rad__mem % alb_nir_can % ptr ( ics : ice , iblk )
      alb_vis_soil => rad__mem % alb_vis_soil % ptr ( ics : ice , iblk )
      alb_nir_soil => rad__mem % alb_nir_soil % ptr ( ics : ice , iblk )
      IF ( use_alb_veg_simple ) THEN
!$acc data &
!$acc present(t_unfilt(:),lai(:),alb_background(:),alb_vis_lnd(:) &
!$acc ,fract_snow_can(:),alb_nir_lnd(:),fract_snow_soil(:),fract_forest(:))
!$acc parallel
!$acc loop gang vector
       DO claw_induction_4 = 1 , size ( alb_vis_lnd , 1 )
        alb_vis_lnd ( claw_induction_4 ) = get_surface_albedo_simple (&
         t_unfilt ( claw_induction_4 ) , fract_forest ( claw_induction_4 ) ,&
         fract_snow_soil ( claw_induction_4 ) , fract_snow_can (&
         claw_induction_4 ) , lai ( claw_induction_4 ) , alb_background (&
         claw_induction_4 ) )
        alb_nir_lnd ( claw_induction_4 ) = alb_vis_lnd ( claw_induction_4 )
       END DO
!$acc end parallel
!$acc end data
      ELSE
       IF ( use_alb_soil_scheme ) THEN
        IF ( tile % lcts ( 1 ) % lib_id == 0 ) THEN
        CALL finish ( trim ("mo_rad_interface" //":update_albedo" ) ,&
        "lctlib not available on tile" )
        END IF
        c_ag_sum_1 => carbon__mem % c_ag_sum_1 % ptr ( ics : ice , iblk )
        CALL calc_soil_albedo ( l_day ( : ) , use_alb_mineralsoil_const ,&
         use_alb_soil_organic_int , use_alb_soil_litter , c_ag_sum_1 ( : ) ,&
         c_bg_sum ( : ) , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
         specificleafarea_c , model % lctlib ( tile % lcts ( 1 ) % lib_id ) %&
         albedolittervis , alb_background ( : ) , alb_vis_mineralsoil ( : ) ,&
         alb_nir_mineralsoil ( : ) , alb_vis_soil ( : ) , alb_nir_soil ( : ) )
       ELSE
        IF ( debug_on ( ) .AND. iblk == 1 ) THEN
        CALL message ( trim ("mo_rad_interface" //":update_albedo" ) ,&
        "Soil albedo scheme not used for veg tiles but taken from bc_land_phy.&
        &nc file" )
        END IF
       END IF
       sky_view_fract => rad__mem % sky_view_fract % ptr ( ics : ice , iblk )
       sky_view_fract_stem => rad__mem % sky_view_fract_stem % ptr ( ics : ice&
        , iblk )
       IF ( tile % lcts ( 1 ) % lib_id == 0 ) THEN
        CALL finish ( trim ("mo_rad_interface" //":update_albedo" ) ,&
        "lctlib not available on tile" )
       END IF
       CALL calc_sky_view_fractions ( l_day ( : ) , lai ( : ) , fract_fpc_max&
        ( : ) , cos_zenith_angle ( : ) , model % lctlib ( tile % lcts ( 1 ) %&
        lib_id ) % stemarea , sky_view_fract ( : ) , sky_view_fract_stem ( : )&
        )
       CALL merge_albedos_of_vegtile ( l_day ( : ) , model % lctlib ( tile %&
        lcts ( 1 ) % lib_id ) % forestflag , sky_view_fract ( : ) ,&
        sky_view_fract_stem ( : ) , albedocanopysnow , fract_snow_soil ( : ) ,&
        fract_snow_can ( : ) , alb_vis_snow ( : ) , alb_nir_snow ( : ) ,&
        alb_vis_soil ( : ) , alb_nir_soil ( : ) , alb_vis_can ( : ) ,&
        alb_nir_can ( : ) , alb_vis_lnd ( : ) , alb_nir_lnd ( : ) )
      END IF
     END IF
    END IF
   END IF
  END IF
  IF ( ( .NOT. tile % is_lake ) ) THEN
!$acc data present(alb_vis(:),alb_vis_lnd(:),alb_nir_lnd(:),alb_nir(:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_5 = 1 , size ( alb_vis , 1 )
    alb_vis ( claw_induction_5 ) = alb_vis_lnd ( claw_induction_5 )
    alb_nir ( claw_induction_5 ) = alb_nir_lnd ( claw_induction_5 )
   END DO
!$acc end parallel
!$acc end data
  END IF
  IF ( trim ( use_alb_soil_organic_c ) =="" ) THEN
   IF ( associated ( c_bg_sum ) ) THEN
!$acc exit data delete(c_bg_sum)
    DEALLOCATE ( c_bg_sum )
   END IF
  END IF
!$acc exit data delete(l_day)
  DEALLOCATE ( l_day )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":update_albedo" ) ,"Finished." )
  END IF
 END SUBROUTINE update_albedo

 SUBROUTINE aggregate_albedo ( tile , options )
  USE mo_jsb_time , ONLY: is_time_experiment_start , is_time_ltrig_rad_m1
  USE mo_jsb_control , ONLY: jsbach_runs_standalone

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  CLASS ( t_jsb_aggregator ) , POINTER :: weighted_by_fract
  TYPE ( t_jsb_model ) , POINTER :: model
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_rad_interface" //&
  ":aggregate_albedo"
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":aggregate_albedo" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  IF ( ( .NOT. is_time_experiment_start ( ) ) .AND. ( .NOT.&
   is_time_ltrig_rad_m1 ( model % id ) ) .AND. ( .NOT. jsbach_runs_standalone&
   ( ) ) ) THEN
   RETURN
  END IF
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  weighted_by_fract => tile % get_aggregator ("weighted_by_fract" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_vis_lnd , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_nir_lnd , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_vis_soil , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_nir_soil , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_vis_mineralsoil ,&
   ics , ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_nir_mineralsoil ,&
   ics , ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_vis_snow , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_nir_snow , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % albedo_lwtr , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % albedo_lice , ics ,&
   ice , iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_vis , ics , ice ,&
   iblk ,"mo_rad_interface" //":aggregate_albedo" )
  CALL weighted_by_fract % aggregate ( tile , rad__mem % alb_nir , ics , ice ,&
   iblk ,"mo_rad_interface" //":aggregate_albedo" )
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_rad_interface" //":aggregate_albedo" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE aggregate_albedo

END MODULE mo_rad_interface

