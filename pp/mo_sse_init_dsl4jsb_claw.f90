MODULE mo_sse_init
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message
 USE mo_jsb_control , ONLY: debug_on
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_grid_class , ONLY: t_jsb_vgrid
 USE mo_jsb_grid , ONLY: get_vgrid
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_io_netcdf , ONLY: t_input_file , jsb_netcdf_open_input
 USE mo_jsb_process_class , ONLY: sse_ , seb_ , hydro_
 USE mo_sse_config_class , ONLY: t_sse_config
 USE mo_hydro_config_class , ONLY: t_hydro_config
 USE mo_sse_memory_class , ONLY: t_sse_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 PRIVATE :: sse_init_bc
 PRIVATE :: sse_init_ic
 PRIVATE :: init_soil_properties
 PUBLIC :: sse_init
 REAL ( KIND= wp ) , PRIVATE , PARAMETER :: fao_thermal_diff ( 0 : 6 ) = (/&
  7.4e-7_wp , 8.7e-7_wp , 8.0e-7_wp , 7.4e-7_wp , 7.1e-7_wp , 6.7e-7_wp , ( -&
  1._wp ) /)
 REAL ( KIND= wp ) , PRIVATE , PARAMETER :: fao_vol_hcap ( 0 : 6 ) = (/&
  2.25e+6_wp , 1.93e+6_wp , 2.10e+6_wp , 2.25e+6_wp , 2.36e+6_wp , 2.48e+6_wp&
  , ( - 1._wp ) /)
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_sse_init"

CONTAINS
 SUBROUTINE sse_init ( tile )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_init" //":sse_init"

  CALL sse_init_bc ( tile )
  CALL sse_init_ic ( tile )
 END SUBROUTINE sse_init

 SUBROUTINE sse_init_bc ( tile )

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  REAL ( KIND= wp ) , POINTER :: fao ( : , : )
  REAL ( KIND= wp ) , POINTER :: ptr_2d ( : , : )
  INTEGER :: iblk
  TYPE ( t_input_file ) :: input_file
  TYPE ( t_jsb_model ) , POINTER :: model
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_init" //":sse_init_bc"

  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( sse_ ) % p % config )
  TYPE IS ( t_sse_config )
   sse__conf => cc
  END SELECT
  IF ( ( .NOT. tile % is_process_active ( sse_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) ) THEN
   CALL message ( trim ("mo_sse_init" //":sse_init_bc" ) ,&
   "Setting  soil energy boundary conditions for tile " // trim ( tile % name&
    ) )
  END IF
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  input_file = jsb_netcdf_open_input ( trim ( sse__conf % bc_filename ) ,&
   model % grid_id )
  fao => input_file % read_2d ( variable_name ="fao" )
  fao ( : , : ) = merge ( fao ( : , : ) , 0._wp , fao ( : , : ) >= 0._wp )
  fao ( : , : ) = max ( real ( lbound ( fao_thermal_diff , 1 ) , wp ) , min (&
   real ( ubound ( fao_thermal_diff , 1 ) , wp ) , fao ( : , : ) ) )
  DO iblk = 1 , size ( fao , 2 ) , 1
   sse__mem % thermal_diffusivity % ptr ( : , iblk ) = fao_thermal_diff ( nint&
    ( fao ( : , iblk ) ) )
  END DO
  IF ( sse__conf % l_heat_cap_map ) THEN
   ptr_2d => input_file % read_2d ( variable_name ="heat_capacity" ,&
    fill_array = sse__mem % vol_heat_cap % ptr )
   ptr_2d ( : , : ) = merge ( ptr_2d ( : , : ) , 2.e6_wp , ptr_2d ( : , : ) >&
    0._wp )
  ELSE
   DO iblk = 1 , size ( fao , 2 ) , 1
    sse__mem % vol_heat_cap % ptr ( : , iblk ) = fao_vol_hcap ( nint ( fao ( :&
     , iblk ) ) )
   END DO
  END IF
  IF ( sse__conf % l_heat_cond_map ) THEN
   ptr_2d => input_file % read_2d ( variable_name ="heat_conductivity" ,&
    fill_array = sse__mem % heat_cond % ptr )
   ptr_2d ( : , : ) = merge ( ptr_2d ( : , : ) , 0.2_wp , ptr_2d ( : , : ) >&
    0._wp )
  ELSE
   sse__mem % heat_cond % ptr ( : , : ) = sse__mem % vol_heat_cap % ptr ( : ,&
    : ) * sse__mem % thermal_diffusivity % ptr ( : , : )
  END IF
!$acc data present(sse__mem%thermal_diffusivity%ptr(:,:), &
!$acc sse__mem%vol_heat_cap%ptr(:,:),  sse__mem%heat_cond%ptr(:,:))
!$acc update device(sse__mem%thermal_diffusivity%ptr(:,:), &
!$acc sse__mem%vol_heat_cap%ptr(:,:),  sse__mem%heat_cond%ptr(:,:))
!$ACC end data
  CALL input_file % close ( )
 END SUBROUTINE sse_init_bc

 SUBROUTINE sse_init_ic ( tile )
  USE mo_sse_process , ONLY: init_soil_temperature

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  REAL ( KIND= wp ) , POINTER :: t_srf_mon ( : , : , : )
  TYPE ( t_input_file ) :: input_file
  TYPE ( t_jsb_vgrid ) , POINTER :: soil_e
  INTEGER :: iblk
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_init" //":sse_init_ic"

  model => get_model ( tile % owner_model_id )
  soil_e => get_vgrid ("soil_depth_energy" )
  SELECT TYPE ( mm => tile % mem ( sse_ ) % p )
  TYPE IS ( t_sse_memory )
   sse__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( cc => model % processes ( sse_ ) % p % config )
  TYPE IS ( t_sse_config )
   sse__conf => cc
  END SELECT
  IF ( ( .NOT. tile % is_process_active ( sse_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) ) THEN
   CALL message ( trim ("mo_sse_init" //":sse_init_ic" ) ,&
   "Initializing soil energy memory for tile " // trim ( tile % name ) )
  END IF
  input_file = jsb_netcdf_open_input ( sse__conf % ic_filename , model %&
   grid_id )
  t_srf_mon => input_file % read_2d_time ( variable_name ="surf_temp" )
  t_srf_mon ( : , : , : ) = merge ( t_srf_mon ( : , : , : ) , 280._wp ,&
   t_srf_mon ( : , : , : ) > 0._wp )
  DO iblk = 1 , size ( t_srf_mon , dim = 2 ) , 1
   CALL init_soil_temperature ( soil_e % levels ( : ) , t_srf_mon ( : , iblk ,&
    1 : 12 ) , sse__mem % t_soil % ptr ( : , : , iblk ) , seb__mem % t % ptr (&
    : , iblk ) )
  END DO
  seb__mem % t_unfilt % ptr ( : , : ) = seb__mem % t % ptr ( : , : )
!$acc data present(seb__mem%t_unfilt%ptr(:,:),  seb__mem%t%ptr(:,:), &
!$acc sse__mem%t_soil%ptr(:,:,:))
!$acc update device(seb__mem%t_unfilt%ptr(:,:),  seb__mem%t%ptr(:,:), &
!$acc sse__mem%t_soil%ptr(:,:,:))
!$ACC end data
  CALL input_file % close ( )
  CALL init_soil_properties ( tile )
 END SUBROUTINE sse_init_ic

 SUBROUTINE init_soil_properties ( tile )
  USE mo_jsb_physical_constants , ONLY: ci , dens_snow_min
  USE mo_sse_constants , ONLY: vol_hcap_ice , vol_hcap_snow , hcond_ice ,&
   hcond_snow
  USE mo_sse_process , ONLY: calc_vol_heat_capacity ,&
   calc_thermal_conductivity

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_sse_config ) , POINTER :: sse__conf
  TYPE ( t_sse_memory ) , POINTER :: sse__mem
  TYPE ( t_hydro_config ) , POINTER :: hydro__conf
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap ( : , : )
  REAL ( KIND= wp ) , POINTER :: heat_cond ( : , : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap_snow ( : , : )
  REAL ( KIND= wp ) , POINTER :: heat_cond_snow ( : , : )
  REAL ( KIND= wp ) , POINTER :: soil_depth_sl ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: fract_org_sl ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap_sl ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: heat_cond_sl ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: w_soil_sl ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: w_ice_sl ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: vol_porosity_sl ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: vol_field_cap_sl ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: fract_water ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: fract_ice ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: soil_e_dz_tmp ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: vol_heat_cap_tmp ( : , : , : )
  REAL ( KIND= wp ) , POINTER :: heat_cond_tmp ( : , : , : )
  TYPE ( t_jsb_model ) , POINTER :: model
  TYPE ( t_jsb_vgrid ) , POINTER :: soil_e
  INTEGER :: nsoil
  INTEGER :: nproma
  INTEGER :: nblks
  INTEGER :: i
  INTEGER :: j
  INTEGER :: b
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_init" //&
  ":init_soil_properties"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4
  INTEGER :: claw_induction_5
  INTEGER :: claw_induction_6
  INTEGER :: claw_induction_7
  INTEGER :: claw_induction_8
  INTEGER :: claw_induction_9
  INTEGER :: claw_induction_10
  INTEGER :: claw_induction_11
  INTEGER :: claw_induction_12
  INTEGER :: claw_induction_13
  INTEGER :: claw_induction_14
  INTEGER :: claw_induction_15
  INTEGER :: claw_induction_16
  INTEGER :: claw_induction_17
  INTEGER :: claw_induction_18
  INTEGER :: claw_induction_19
  INTEGER :: claw_induction_20
  INTEGER :: claw_induction_21
  INTEGER :: claw_induction_22
  INTEGER :: claw_induction_23
  INTEGER :: claw_induction_24
  INTEGER :: claw_induction_25
  INTEGER :: claw_induction_26
  INTEGER :: claw_induction_27
  INTEGER :: claw_induction_28
  INTEGER :: claw_induction_29
  INTEGER :: claw_induction_30

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
  vol_heat_cap_sl => sse__mem % vol_heat_cap_sl % ptr ( : , : , : )
  heat_cond_sl => sse__mem % heat_cond_sl % ptr ( : , : , : )
  vol_heat_cap_snow => sse__mem % vol_heat_cap_snow % ptr ( : , : )
  heat_cond_snow => sse__mem % heat_cond_snow % ptr ( : , : )
  IF ( tile % is_glacier ) THEN
!$acc data present(vol_heat_cap_sl(:,:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_16 = 1 , size ( vol_heat_cap_sl , 1 )
    DO claw_induction_17 = 1 , size ( vol_heat_cap_sl , 2 )
     DO claw_induction_18 = 1 , size ( vol_heat_cap_sl , 3 )
      vol_heat_cap_sl ( claw_induction_16 , claw_induction_17 ,&
       claw_induction_18 ) = vol_hcap_ice
     END DO
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc data present(heat_cond_sl(:,:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_13 = 1 , size ( heat_cond_sl , 1 )
    DO claw_induction_14 = 1 , size ( heat_cond_sl , 2 )
     DO claw_induction_15 = 1 , size ( heat_cond_sl , 3 )
      heat_cond_sl ( claw_induction_13 , claw_induction_14 , claw_induction_15&
       ) = hcond_ice
     END DO
    END DO
   END DO
!$acc end parallel
!$acc end data
!$cla$w end expand
!$acc data present(vol_heat_cap_snow(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_11 = 1 , size ( vol_heat_cap_snow , 1 )
    DO claw_induction_12 = 1 , size ( vol_heat_cap_snow , 2 )
     vol_heat_cap_snow ( claw_induction_11 , claw_induction_12 ) =&
      vol_hcap_snow
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc data present(heat_cond_snow(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_9 = 1 , size ( heat_cond_snow , 1 )
    DO claw_induction_10 = 1 , size ( heat_cond_snow , 2 )
     heat_cond_snow ( claw_induction_9 , claw_induction_10 ) = hcond_snow
    END DO
   END DO
!$acc end parallel
!$acc end data
!$cla$w end expand
   RETURN
  END IF
  w_soil_sl => hydro__mem % w_soil_sl % ptr ( : , : , : )
  w_ice_sl => hydro__mem % w_ice_sl % ptr ( : , : , : )
  vol_porosity_sl => hydro__mem % vol_porosity_sl % ptr ( : , : , : )
  vol_field_cap_sl => hydro__mem % vol_field_cap_sl % ptr ( : , : , : )
  soil_depth_sl => hydro__mem % soil_depth_sl % ptr ( : , : , : )
  vol_heat_cap => sse__mem % vol_heat_cap % ptr ( : , : )
  heat_cond => sse__mem % heat_cond % ptr ( : , : )
  nproma = size ( heat_cond , 1 )
  nblks = size ( heat_cond , 2 )
  IF ( hydro__conf % l_organic .AND. tile % contains_soil ) THEN
   fract_org_sl => hydro__mem % fract_org_sl % ptr ( : , : , : )
  ELSE
   ALLOCATE ( fract_org_sl ( nproma , nsoil , nblks ) )
!$acc enter data create(fract_org_sl)
!$acc data present(fract_org_sl(:,:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_6 = 1 , size ( fract_org_sl , 1 )
    DO claw_induction_7 = 1 , size ( fract_org_sl , 2 )
     DO claw_induction_8 = 1 , size ( fract_org_sl , 3 )
      fract_org_sl ( claw_induction_6 , claw_induction_7 , claw_induction_8 )&
       = 0._wp
     END DO
    END DO
   END DO
!$acc end parallel
!$acc end data
  END IF
  ALLOCATE ( soil_e_dz_tmp ( nproma , nsoil , nblks ) , vol_heat_cap_tmp (&
   nproma , nsoil , nblks ) , heat_cond_tmp ( nproma , nsoil , nblks ) )
!$acc enter data pcreate(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)
  soil_e_dz_tmp ( : , : , : ) = spread ( spread ( soil_e % dz ( : ) , dim = 1&
   , ncopies = nproma ) , dim = 3 , ncopies = nblks )
  vol_heat_cap_tmp ( : , : , : ) = spread ( vol_heat_cap ( : , : ) , dim = 2 ,&
   ncopies = nsoil )
  heat_cond_tmp ( : , : , : ) = spread ( heat_cond ( : , : ) , dim = 2 ,&
   ncopies = nsoil )
!$acc update device(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)
  IF ( sse__conf % l_heat_cap_dyn ) THEN
   ALLOCATE ( fract_water ( nproma , nsoil , nblks ) , fract_ice ( nproma ,&
    nsoil , nblks ) )
!$acc enter data pcreate(fract_water,fract_ice)
!$acc data present(soil_depth_sl,fract_water,fract_ice,w_soil_sl,w_ice_sl)
!$acc update device(soil_depth_sl,w_soil_sl,w_ice_sl)
!$acc parallel loop gang vector collapse(3)
   DO i = 1 , nproma , 1
    DO j = 1 , nsoil , 1
     DO b = 1 , nblks , 1
      IF ( soil_depth_sl ( i , j , b ) > 0._wp ) THEN
       fract_water ( i , j , b ) = w_soil_sl ( i , j , b ) / soil_depth_sl ( i&
        , j , b )
       fract_ice ( i , j , b ) = w_ice_sl ( i , j , b ) / soil_depth_sl ( i ,&
        j , b )
      ELSE
       fract_water ( i , j , b ) = 0._wp
       fract_ice ( i , j , b ) = 0._wp
      END IF
     END DO
    END DO
   END DO
!$acc update host(fract_water,fract_ice)
!$acc end data
!$acc data &
!$acc present(vol_porosity_sl(:,:,:),soil_e_dz_tmp(:,:,:),fract_ice(:,:,:) &
!$acc ,fract_water(:,:,:),soil_depth_sl(:,:,:),vol_heat_cap_tmp(:,:,:) &
!$acc ,vol_heat_cap_sl(:,:,:),fract_org_sl(:,:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_3 = 1 , size ( vol_heat_cap_sl , 1 )
    DO claw_induction_4 = 1 , size ( vol_heat_cap_sl , 2 )
     DO claw_induction_5 = 1 , size ( vol_heat_cap_sl , 3 )
      vol_heat_cap_sl ( claw_induction_3 , claw_induction_4 , claw_induction_5&
       ) = calc_vol_heat_capacity ( vol_heat_cap_tmp ( claw_induction_3 ,&
       claw_induction_4 , claw_induction_5 ) , vol_porosity_sl (&
       claw_induction_3 , claw_induction_4 , claw_induction_5 ) , fract_org_sl&
       ( claw_induction_3 , claw_induction_4 , claw_induction_5 ) ,&
       fract_water ( claw_induction_3 , claw_induction_4 , claw_induction_5 )&
       , fract_ice ( claw_induction_3 , claw_induction_4 , claw_induction_5 )&
       , soil_depth_sl ( claw_induction_3 , claw_induction_4 ,&
       claw_induction_5 ) , soil_e_dz_tmp ( claw_induction_3 ,&
       claw_induction_4 , claw_induction_5 ) )
     END DO
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc exit data delete(fract_water, fract_ice)
   DEALLOCATE ( fract_water , fract_ice )
  ELSE
!$acc data present(vol_heat_cap_sl,vol_heat_cap)
!$acc update device(vol_heat_cap)
!$acc parallel loop gang vector collapse(3)
   DO i = 1 , nproma , 1
    DO j = 1 , nsoil , 1
     DO b = 1 , nblks , 1
      vol_heat_cap_sl ( i , j , b ) = vol_heat_cap ( i , b )
     END DO
    END DO
   END DO
!$acc update host(vol_heat_cap_sl)
!$acc end data
  END IF
  IF ( sse__conf % l_heat_cond_dyn ) THEN
!$acc data &
!$acc present(w_soil_sl(:,:,:),vol_porosity_sl(:,:,:),soil_e_dz_tmp(:,:,:) &
!$acc ,heat_cond_tmp(:,:,:),w_ice_sl(:,:,:),soil_depth_sl(:,:,:) &
!$acc ,heat_cond_sl(:,:,:),vol_field_cap_sl(:,:,:),fract_org_sl(:,:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_0 = 1 , size ( heat_cond_sl , 1 )
    DO claw_induction_1 = 1 , size ( heat_cond_sl , 2 )
     DO claw_induction_2 = 1 , size ( heat_cond_sl , 3 )
      heat_cond_sl ( claw_induction_0 , claw_induction_1 , claw_induction_2 )&
       = calc_thermal_conductivity ( heat_cond_tmp ( claw_induction_0 ,&
       claw_induction_1 , claw_induction_2 ) , w_soil_sl ( claw_induction_0 ,&
       claw_induction_1 , claw_induction_2 ) , w_ice_sl ( claw_induction_0 ,&
       claw_induction_1 , claw_induction_2 ) , vol_porosity_sl (&
       claw_induction_0 , claw_induction_1 , claw_induction_2 ) ,&
       vol_field_cap_sl ( claw_induction_0 , claw_induction_1 ,&
       claw_induction_2 ) , fract_org_sl ( claw_induction_0 , claw_induction_1&
       , claw_induction_2 ) , soil_depth_sl ( claw_induction_0 ,&
       claw_induction_1 , claw_induction_2 ) , soil_e_dz_tmp (&
       claw_induction_0 , claw_induction_1 , claw_induction_2 ) )
     END DO
    END DO
   END DO
!$acc end parallel
!$acc end data
  ELSE
!$acc data present(heat_cond_sl,heat_cond)
!$acc update device(heat_cond)
!$acc parallel loop gang vector collapse(3)
   DO i = 1 , nproma , 1
    DO j = 1 , nsoil , 1
     DO b = 1 , nblks , 1
      heat_cond_sl ( i , j , b ) = heat_cond ( i , b )
     END DO
    END DO
   END DO
!$acc update host(heat_cond_sl)
!$acc end data
  END IF
  IF ( ( .NOT. ( hydro__conf % l_organic .AND. tile % contains_soil ) ) ) THEN
!$acc exit data delete(fract_org_sl)
   DEALLOCATE ( fract_org_sl )
  END IF
  IF ( sse__conf % l_snow ) THEN
   IF ( model % config % l_compat401 ) THEN
!$acc data present(vol_heat_cap_snow(:,:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_29 = 1 , size ( vol_heat_cap_snow , 1 )
     DO claw_induction_30 = 1 , size ( vol_heat_cap_snow , 2 )
      vol_heat_cap_snow ( claw_induction_29 , claw_induction_30 ) = 634500._wp
     END DO
    END DO
!$acc end parallel
!$acc end data
!$acc data present(heat_cond_snow(:,:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_27 = 1 , size ( heat_cond_snow , 1 )
     DO claw_induction_28 = 1 , size ( heat_cond_snow , 2 )
      heat_cond_snow ( claw_induction_27 , claw_induction_28 ) = 0.31_wp
     END DO
    END DO
!$acc end parallel
!$acc end data
!$cla$w end expand
   ELSE
    IF ( sse__conf % l_dynsnow ) THEN
!$acc data present(vol_heat_cap_snow(:,:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_25 = 1 , size ( vol_heat_cap_snow , 1 )
      DO claw_induction_26 = 1 , size ( vol_heat_cap_snow , 2 )
       vol_heat_cap_snow ( claw_induction_25 , claw_induction_26 ) = ci *&
        dens_snow_min
      END DO
     END DO
!$acc end parallel
!$acc end data
!$acc data present(heat_cond_snow(:,:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_23 = 1 , size ( heat_cond_snow , 1 )
      DO claw_induction_24 = 1 , size ( heat_cond_snow , 2 )
       heat_cond_snow ( claw_induction_23 , claw_induction_24 ) = 2.5e-6_wp *&
        dens_snow_min ** ( 2._wp ) - 0.123e-3_wp * dens_snow_min + 0.024_wp
      END DO
     END DO
!$acc end parallel
!$acc end data
!$cla$w end expand
    ELSE
!$acc data present(vol_heat_cap_snow(:,:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_21 = 1 , size ( vol_heat_cap_snow , 1 )
      DO claw_induction_22 = 1 , size ( vol_heat_cap_snow , 2 )
       vol_heat_cap_snow ( claw_induction_21 , claw_induction_22 ) =&
        vol_hcap_snow
      END DO
     END DO
!$acc end parallel
!$acc end data
!$acc data present(heat_cond_snow(:,:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_19 = 1 , size ( heat_cond_snow , 1 )
      DO claw_induction_20 = 1 , size ( heat_cond_snow , 2 )
       heat_cond_snow ( claw_induction_19 , claw_induction_20 ) = hcond_snow
      END DO
     END DO
!$acc end parallel
!$acc end data
!$cla$w end expand
    END IF
   END IF
  END IF
!$ACC exit data delete(soil_e_dz_tmp,vol_heat_cap_tmp,heat_cond_tmp)
  DEALLOCATE ( soil_e_dz_tmp , vol_heat_cap_tmp , heat_cond_tmp )
 END SUBROUTINE init_soil_properties

END MODULE mo_sse_init

