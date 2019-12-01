MODULE mo_pheno_init
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message
 USE mo_jsb_control , ONLY: debug_on
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract
 USE mo_jsb_io_netcdf , ONLY: t_input_file , jsb_netcdf_open_input
 USE mo_jsb_process_class , ONLY: pheno_
 USE mo_pheno_config_class , ONLY: t_pheno_config
 USE mo_pheno_memory_class , ONLY: t_pheno_memory
 PRIVATE :: pheno_init_bc
 PRIVATE :: pheno_init_ic
 PUBLIC :: pheno_init
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_pheno_init"

CONTAINS
 SUBROUTINE pheno_init ( tile )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_init" //":pheno_init"

  CALL pheno_init_bc ( tile )
  CALL pheno_init_ic ( tile )
 END SUBROUTINE pheno_init

 SUBROUTINE pheno_init_bc ( tile )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_pheno_config ) , POINTER :: pheno__conf
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  REAL ( KIND= wp ) , POINTER :: ptr_mon ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: ptr_2d ( : , : )
  TYPE ( t_input_file ) :: input_file
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_init" //&
  ":pheno_init_bc"

  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( pheno_ ) % p % config )
  TYPE IS ( t_pheno_config )
   pheno__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  IF ( debug_on ( ) ) THEN
   CALL message ( trim ("mo_pheno_init" //":pheno_init_bc" ) ,&
   "Setting boundary conditions of pheno memory for tile " // trim ( tile %&
    name ) //" from " // trim ( pheno__conf % bc_filename ) )
  END IF
  input_file = jsb_netcdf_open_input ( trim ( pheno__conf % bc_filename ) ,&
   model % grid_id )
  IF ( debug_on ( ) ) THEN
   CALL message ( trim ("mo_pheno_init" //":pheno_init_bc" ) ,&
   "reading lai_clim ..." )
  END IF
  ptr_mon => input_file % read_2d_time ( variable_name ="lai_clim" )
  ptr_mon ( : , : , : ) = merge ( ptr_mon ( : , : , : ) , 0._wp , ptr_mon ( :&
   , : , : ) > 0._wp )
  pheno__mem % lai_mon ( : , : , 1 : 12 ) = ptr_mon ( : , : , 1 : 12 )
  pheno__mem % lai_mon ( : , : , 0 ) = ptr_mon ( : , : , 12 )
  pheno__mem % lai_mon ( : , : , 13 ) = ptr_mon ( : , : , 1 )
  DEALLOCATE ( ptr_mon )
!$ACC enter data attach(pheno__mem%lai_mon)
!$ACC update device(pheno__mem%lai_mon)
  IF ( debug_on ( ) ) THEN
   CALL message ( trim ("mo_pheno_init" //":pheno_init_bc" ) ,&
   "reading veg_fract ..." )
  END IF
  ptr_mon => input_file % read_2d_time ( variable_name ="veg_fract" )
  ptr_mon ( : , : , : ) = merge ( ptr_mon ( : , : , : ) , 0._wp , ptr_mon ( :&
   , : , : ) > 0._wp )
  pheno__mem % fract_fpc_mon ( : , : , 1 : 12 ) = ptr_mon ( : , : , 1 : 12 )
  pheno__mem % fract_fpc_mon ( : , : , 0 ) = ptr_mon ( : , : , 12 )
  pheno__mem % fract_fpc_mon ( : , : , 13 ) = ptr_mon ( : , : , 1 )
  DEALLOCATE ( ptr_mon )
!$ACC update device(pheno__mem%fract_fpc_mon)
  IF ( debug_on ( ) ) THEN
   CALL message ( trim ("mo_pheno_init" //":pheno_init_bc" ) ,&
   "reading forest_fract ..." )
  END IF
  ptr_2d => input_file % read_2d ( variable_name ="forest_fract" , fill_array&
   = pheno__mem % fract_forest % ptr )
  pheno__mem % fract_forest % ptr ( : , : ) = merge ( pheno__mem %&
   fract_forest % ptr ( : , : ) , 0._wp , pheno__mem % fract_forest % ptr ( :&
   , : ) >= 0._wp )
!$ACC update device(pheno__mem%fract_forest%ptr(:,:))
  CALL input_file % close ( )
  input_file = jsb_netcdf_open_input ( trim ( model % config % fract_filename&
   ) , model % grid_id )
  IF ( debug_on ( ) ) THEN
   CALL message ( trim ("mo_pheno_init" //":pheno_init_bc" ) ,&
   "reading veg_ratio_max ..." )
  END IF
  ptr_2d => input_file % read_2d ( variable_name ="veg_ratio_max" , fill_array&
   = pheno__mem % fract_fpc_max % ptr )
  pheno__mem % fract_fpc_max % ptr ( : , : ) = merge ( pheno__mem %&
   fract_fpc_max % ptr ( : , : ) , 0._wp , pheno__mem % fract_fpc_max % ptr (&
   : , : ) >= 0._wp )
!$ACC update device(pheno__mem%fract_fpc_max%ptr(:,:))
  CALL input_file % close ( )
  IF ( associated ( ptr_2d ) ) THEN
   CONTINUE
  END IF
  NULLIFY ( ptr_2d )
  IF ( tile % lcts ( 1 ) % lib_id /= 0 ) THEN
   IF ( pheno__conf % l_forestregrowth .AND. ( .NOT. model % lctlib ( tile %&
    lcts ( 1 ) % lib_id ) % forestflag ) ) THEN
    pheno__mem % maxlai_allom % ptr ( : , : ) = model % lctlib ( tile % lcts (&
     1 ) % lib_id ) % maxlai
!$ACC update device(pheno__mem%maxlai_allom%ptr(:,:))
   END IF
   IF ( pheno__conf % l_forestregrowth .AND. model % lctlib ( tile % lcts ( 1&
    ) % lib_id ) % forestflag ) THEN
    pheno__mem % veg_fract_correction % ptr ( : , : ) = 1.0_wp - exp ( ( -&
     pheno__mem % maxlai_allom % ptr ( : , : ) / model % lctlib ( tile % lcts&
     ( 1 ) % lib_id ) % clumpinessfactor ) )
   ELSE
    pheno__mem % veg_fract_correction % ptr ( : , : ) = 1.0_wp - exp ( ( -&
     model % lctlib ( tile % lcts ( 1 ) % lib_id ) % maxlai / model % lctlib (&
     tile % lcts ( 1 ) % lib_id ) % clumpinessfactor ) )
   END IF
!$ACC update device(pheno__mem%veg_fract_correction%ptr(:,:))
  END IF
 END SUBROUTINE pheno_init_bc

 SUBROUTINE pheno_init_ic ( tile )
  USE mo_pheno_process , ONLY: get_foliage_projected_cover
  USE mo_jsb_time , ONLY: get_time_interpolation_weights

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_pheno_config ) , POINTER :: pheno__conf
  TYPE ( t_pheno_memory ) , POINTER :: pheno__mem
  REAL ( KIND= wp ) :: wgt1
  REAL ( KIND= wp ) :: wgt2
  INTEGER :: nmw1
  INTEGER :: nmw2
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_pheno_init" //&
  ":pheno_init_ic"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_3

  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( pheno_ ) % p % config )
  TYPE IS ( t_pheno_config )
   pheno__conf => cc
  END SELECT
  SELECT TYPE ( mm => tile % mem ( pheno_ ) % p )
  TYPE IS ( t_pheno_memory )
   pheno__mem => mm
  END SELECT
  IF ( debug_on ( ) ) THEN
   CALL message ( trim ("mo_pheno_init" //":pheno_init_ic" ) ,&
   "Setting initial conditions of pheno memory (lai) for tile " // trim ( tile&
    % name ) //" from " // trim ( pheno__conf % bc_filename ) )
  END IF
  CALL get_time_interpolation_weights ( wgt1 , wgt2 , nmw1 , nmw2 )
  pheno__mem % lai % ptr ( : , : ) = wgt1 * pheno__mem % lai_mon ( : , : ,&
   nmw1 ) + wgt2 * pheno__mem % lai_mon ( : , : , nmw2 )
  IF ( model % config % l_compat401 ) THEN
   DO claw_induction_0 = 1 , size ( pheno__mem % fract_fpc % ptr , 1 )
    DO claw_induction_1 = 1 , size ( pheno__mem % fract_fpc % ptr , 2 )
     pheno__mem % fract_fpc % ptr ( claw_induction_0 , claw_induction_1 ) =&
      get_foliage_projected_cover ( pheno__mem % fract_fpc_max % ptr (&
      claw_induction_0 , claw_induction_1 ) , pheno__mem % lai % ptr (&
      claw_induction_0 , claw_induction_1 ) , 2._wp )
    END DO
   END DO
!$ACC update device(pheno__mem%fract_fpc%ptr(:,:))
  ELSE
   IF ( trim ( pheno__conf % scheme ) =="climatology" .OR. tile % lcts ( 1 ) %&
    lib_id == 0 ) THEN
    pheno__mem % fract_fpc % ptr ( : , : ) = wgt1 * pheno__mem % fract_fpc_mon&
     ( : , : , nmw1 ) + wgt2 * pheno__mem % fract_fpc_mon ( : , : , nmw2 )
   ELSE
    IF ( ( .NOT. pheno__conf % l_forestregrowth ) ) THEN
     pheno__mem % lai % ptr ( : , : ) = min ( pheno__mem % lai % ptr ( : , : )&
      , model % lctlib ( tile % lcts ( 1 ) % lib_id ) % maxlai )
    ELSE
     pheno__mem % lai % ptr ( : , : ) = min ( pheno__mem % lai % ptr ( : , : )&
      , pheno__mem % maxlai_allom % ptr ( : , : ) )
    END IF
    DO claw_induction_2 = 1 , size ( pheno__mem % fract_fpc % ptr , 1 )
     DO claw_induction_3 = 1 , size ( pheno__mem % fract_fpc % ptr , 2 )
      pheno__mem % fract_fpc % ptr ( claw_induction_2 , claw_induction_3 ) =&
       get_foliage_projected_cover ( pheno__mem % fract_fpc_max % ptr (&
       claw_induction_2 , claw_induction_3 ) , pheno__mem % lai % ptr (&
       claw_induction_2 , claw_induction_3 ) , model % lctlib ( tile % lcts (&
       1 ) % lib_id ) % clumpinessfactor )
     END DO
    END DO
   END IF
!$ACC update device(pheno__mem%fract_fpc%ptr(:,:))
  END IF
!$ACC update device(pheno__mem%lai%ptr(:,:))
 END SUBROUTINE pheno_init_ic

END MODULE mo_pheno_init

