MODULE mo_seb_lake
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: message , finish
 USE mo_jsb_model_class , ONLY: t_jsb_model
 USE mo_jsb_class , ONLY: get_model
 USE mo_jsb_tile_class , ONLY: t_jsb_tile_abstract
 USE mo_jsb_task_class , ONLY: t_jsb_task_options
 USE mo_jsb_control , ONLY: debug_on
 USE mo_jsb_process_class , ONLY: seb_ , rad_ , hydro_ , a2l_
 USE mo_seb_config_class , ONLY: t_seb_config
 USE mo_a2l_memory_class , ONLY: t_a2l_memory
 USE mo_seb_memory_class , ONLY: t_seb_memory
 USE mo_rad_memory_class , ONLY: t_rad_memory
 USE mo_hydro_memory_class , ONLY: t_hydro_memory
 PRIVATE :: calc_lake_wtr_temp
 PRIVATE :: calc_lake_ice_temp
 PRIVATE :: calc_lice_thickness
 PUBLIC :: update_surface_energy_lake
 PUBLIC :: update_surface_fluxes_lake
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_seb_lake"

CONTAINS
 SUBROUTINE update_surface_energy_lake ( tile , options )
  USE mo_phy_schemes , ONLY: qsat , surface_dry_static_energy

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_config ) , POINTER :: seb__conf
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  TYPE ( t_rad_memory ) , POINTER :: rad__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  REAL ( KIND= wp ) , POINTER :: t ( : )
  REAL ( KIND= wp ) , POINTER :: t_unfilt ( : )
  REAL ( KIND= wp ) , POINTER :: t_lwtr ( : )
  REAL ( KIND= wp ) , POINTER :: t_lice ( : )
  REAL ( KIND= wp ) , POINTER :: qsat_lwtr ( : )
  REAL ( KIND= wp ) , POINTER :: qsat_lice ( : )
  REAL ( KIND= wp ) , POINTER :: s_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: s_ice ( : )
  REAL ( KIND= wp ) , POINTER :: heat_cap ( : )
  REAL ( KIND= wp ) , POINTER :: press_srf ( : )
  REAL ( KIND= wp ) , POINTER :: snow ( : )
  REAL ( KIND= wp ) , POINTER :: latent_hflx_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: sensible_hflx_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: latent_hflx_ice ( : )
  REAL ( KIND= wp ) , POINTER :: sensible_hflx_ice ( : )
  REAL ( KIND= wp ) , POINTER :: hflx_form_ice ( : )
  REAL ( KIND= wp ) , POINTER :: hflx_cond_ice ( : )
  REAL ( KIND= wp ) , POINTER :: hflx_melt_lice ( : )
  REAL ( KIND= wp ) , POINTER :: hflx_melt_snow ( : )
  REAL ( KIND= wp ) , POINTER :: fract_lice ( : )
  REAL ( KIND= wp ) , POINTER :: depth_lice ( : )
  REAL ( KIND= wp ) , POINTER :: fract_snow_lice ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow_lice ( : )
  REAL ( KIND= wp ) , POINTER :: w_snow ( : )
  REAL ( KIND= wp ) , POINTER :: rad_net_lwtr ( : )
  REAL ( KIND= wp ) , POINTER :: rad_net_lice ( : )
  REAL ( KIND= wp ) :: hcap_wtr ( 1 : options % nc )
  REAL ( KIND= wp ) :: hcap_ice ( 1 : options % nc )
  REAL ( KIND= wp ) :: zero ( 1 : options % nc )
  REAL ( KIND= wp ) :: lhflx_tmp ( 1 : options % nc )
  REAL ( KIND= wp ) :: hflx_tmp ( 1 : options % nc )
  LOGICAL :: l_ice ( 1 : options % nc )
  LOGICAL :: l_ice_formation ( 1 : options % nc )
  LOGICAL :: l_ice_melting ( 1 : options % nc )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: dtime
  REAL ( KIND= wp ) :: ml_depth
  REAL ( KIND= wp ) :: min_ice_depth
  TYPE ( t_jsb_model ) , POINTER :: model
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_lake" //&
  ":update_surface_energy_lake"
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_3

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  dtime = options % dtime
  IF ( ( .NOT. tile % is_process_active ( seb_ ) ) .OR. ( .NOT. tile %&
   contains_lake ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_lake" //":update_surface_energy_lake" ) ,&
   "Starting on tile " // trim ( tile % name ) //" ..." )
  END IF
  model => get_model ( tile % owner_model_id )
  SELECT TYPE ( cc => model % processes ( seb_ ) % p % config )
  TYPE IS ( t_seb_config )
   seb__conf => cc
  END SELECT
  ml_depth = seb__conf % lake_mixed_layer_depth
  min_ice_depth = seb__conf % lake_min_ice_depth
  SELECT TYPE ( mm => tile % mem ( seb_ ) % p )
  TYPE IS ( t_seb_memory )
   seb__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( 1 ) % p )
  TYPE IS ( t_a2l_memory )
   a2l__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( rad_ ) % p )
  TYPE IS ( t_rad_memory )
   rad__mem => mm
  END SELECT
  SELECT TYPE ( mm => tile % mem ( hydro_ ) % p )
  TYPE IS ( t_hydro_memory )
   hydro__mem => mm
  END SELECT
  t => seb__mem % t % ptr ( ics : ice , iblk )
  t_unfilt => seb__mem % t_unfilt % ptr ( ics : ice , iblk )
  t_lwtr => seb__mem % t_lwtr % ptr ( ics : ice , iblk )
  t_lice => seb__mem % t_lice % ptr ( ics : ice , iblk )
  qsat_lwtr => seb__mem % qsat_lwtr % ptr ( ics : ice , iblk )
  qsat_lice => seb__mem % qsat_lice % ptr ( ics : ice , iblk )
  s_wtr => seb__mem % s_wtr % ptr ( ics : ice , iblk )
  s_ice => seb__mem % s_ice % ptr ( ics : ice , iblk )
  heat_cap => seb__mem % heat_cap % ptr ( ics : ice , iblk )
  hflx_form_ice => seb__mem % hflx_form_ice % ptr ( ics : ice , iblk )
  hflx_cond_ice => seb__mem % hflx_cond_ice % ptr ( ics : ice , iblk )
  hflx_melt_lice => seb__mem % hflx_melt_lice % ptr ( ics : ice , iblk )
  hflx_melt_snow => seb__mem % hflx_melt_snow % ptr ( ics : ice , iblk )
  press_srf => a2l__mem % press_srf % ptr ( ics : ice , iblk )
  snow => a2l__mem % snow % ptr ( ics : ice , iblk )
  latent_hflx_wtr => seb__mem % latent_hflx_wtr % ptr ( ics : ice , iblk )
  sensible_hflx_wtr => seb__mem % sensible_hflx_wtr % ptr ( ics : ice , iblk )
  latent_hflx_ice => seb__mem % latent_hflx_ice % ptr ( ics : ice , iblk )
  sensible_hflx_ice => seb__mem % sensible_hflx_ice % ptr ( ics : ice , iblk )
  fract_lice => seb__mem % fract_lice % ptr ( ics : ice , iblk )
  depth_lice => seb__mem % depth_lice % ptr ( ics : ice , iblk )
  fract_snow_lice => hydro__mem % fract_snow_lice % ptr ( ics : ice , iblk )
  w_snow_lice => hydro__mem % w_snow_lice % ptr ( ics : ice , iblk )
  w_snow => hydro__mem % w_snow % ptr ( ics : ice , iblk )
  rad_net_lwtr => rad__mem % rad_net_lwtr % ptr ( ics : ice , iblk )
  rad_net_lice => rad__mem % rad_net_lice % ptr ( ics : ice , iblk )
!$acc data create(hcap_wtr,hcap_ice,zero,l_ice, &
!$acc l_ice_formation,l_ice_melting,lhflx_tmp,hflx_tmp)
!$acc data &
!$acc present(rad_net_lwtr(:),l_ice_formation(:),l_ice_melting(:),l_ice(:) &
!$acc ,zero(:),hflx_tmp(:),sensible_hflx_wtr(:),latent_hflx_wtr(:) &
!$acc ,hcap_ice(:),fract_lice(:),hcap_wtr(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_0 = 1 , size ( zero , 1 )
   zero ( claw_induction_0 ) = 0._wp
   hcap_wtr ( claw_induction_0 ) = 0._wp
   hcap_ice ( claw_induction_0 ) = 0._wp
   l_ice_formation ( claw_induction_0 ) = .FALSE.
   l_ice_melting ( claw_induction_0 ) = .FALSE.
   l_ice ( claw_induction_0 ) = fract_lice ( claw_induction_0 ) > 0.5_wp
   hflx_tmp ( claw_induction_0 ) = rad_net_lwtr ( claw_induction_0 ) +&
    latent_hflx_wtr ( claw_induction_0 ) + sensible_hflx_wtr (&
    claw_induction_0 )
  END DO
!$acc end parallel
!$acc end data
  CALL calc_lake_wtr_temp ( t_lwtr ( : ) , hflx_tmp ( : ) , hflx_form_ice ( :&
   ) , hcap_wtr ( : ) , l_ice_formation ( : ) , l_ice_melting ( : ) , l_ice (&
   : ) , dtime , ml_depth , min_ice_depth )
  CALL calc_lice_thickness ( depth_lice ( : ) , fract_lice ( : ) ,&
   hflx_form_ice ( : ) , latent_hflx_ice ( : ) , hflx_cond_ice ( : ) ,&
   hflx_melt_lice ( : ) , l_ice_formation ( : ) , l_ice_melting ( : ) ,&
   fract_snow_lice ( : ) , min_ice_depth , dtime )
  CALL calc_lake_wtr_temp ( t_lwtr ( : ) , zero ( : ) , hflx_melt_lice ( : ) ,&
   hcap_wtr ( : ) , l_ice_formation ( : ) , l_ice_melting ( : ) , l_ice ( : )&
   , dtime , ml_depth , min_ice_depth )
!$acc data &
!$acc present(latent_hflx_ice(:),rad_net_lice(:),lhflx_tmp(:) &
!$acc ,sensible_hflx_ice(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_3 = 1 , size ( lhflx_tmp , 1 )
   lhflx_tmp ( claw_induction_3 ) = rad_net_lice ( claw_induction_3 ) +&
    latent_hflx_ice ( claw_induction_3 ) + sensible_hflx_ice (&
    claw_induction_3 )
  END DO
!$acc end parallel
!$acc end data
  CALL calc_lake_ice_temp ( t_lice ( : ) , hcap_ice ( : ) , hflx_cond_ice ( :&
   ) , hflx_melt_snow ( : ) , hflx_melt_lice ( : ) , lhflx_tmp ( : ) ,&
   latent_hflx_ice ( : ) , snow ( : ) , depth_lice ( : ) , w_snow_lice ( : ) ,&
   fract_snow_lice ( : ) , min_ice_depth , dtime )
!$acc data &
!$acc present(t_lice(:),hcap_wtr(:),qsat_lwtr(:),s_wtr(:),t_unfilt(:) &
!$acc ,w_snow_lice(:),qsat_lice(:),t(:),heat_cap(:),t_lwtr(:),press_srf(:) &
!$acc ,s_ice(:),hcap_ice(:),fract_lice(:),w_snow(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_1 = 1 , size ( qsat_lwtr , 1 )
   qsat_lwtr ( claw_induction_1 ) = qsat ( t_lwtr ( claw_induction_1 ) ,&
    press_srf ( claw_induction_1 ) )
   qsat_lice ( claw_induction_1 ) = qsat ( t_lice ( claw_induction_1 ) ,&
    press_srf ( claw_induction_1 ) )
   s_wtr ( claw_induction_1 ) = surface_dry_static_energy ( t_lwtr (&
    claw_induction_1 ) , qsat_lwtr ( claw_induction_1 ) )
   s_ice ( claw_induction_1 ) = surface_dry_static_energy ( t_lice (&
    claw_induction_1 ) , qsat_lice ( claw_induction_1 ) )
   t ( claw_induction_1 ) = ( 1._wp - fract_lice ( claw_induction_1 ) ) *&
    t_lwtr ( claw_induction_1 ) + fract_lice ( claw_induction_1 ) * t_lice (&
    claw_induction_1 )
   t_unfilt ( claw_induction_1 ) = t ( claw_induction_1 )
   heat_cap ( claw_induction_1 ) = ( 1._wp - fract_lice ( claw_induction_1 ) )&
    * hcap_wtr ( claw_induction_1 ) + fract_lice ( claw_induction_1 ) *&
    hcap_ice ( claw_induction_1 )
   w_snow ( claw_induction_1 ) = fract_lice ( claw_induction_1 ) * w_snow_lice&
    ( claw_induction_1 )
  END DO
!$acc end parallel
!$acc end data
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_lake" //":update_surface_energy_lake" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_surface_energy_lake

 SUBROUTINE calc_lake_wtr_temp ( t_lwtr , hflx , hflx_resid , hcap ,&
  l_ice_formation , l_ice_melting , l_ice , dtime , ml_depth , min_ice_depth )
  USE mo_jsb_physical_constants , ONLY: tmelt , rhoh2o , clw , rhoi , alf

  REAL ( KIND= wp ) , INTENT(INOUT) :: t_lwtr ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: hflx ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: hflx_resid ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: hcap ( : )
  LOGICAL , INTENT(INOUT) :: l_ice_formation ( : )
  LOGICAL , INTENT(INOUT) :: l_ice_melting ( : )
  LOGICAL , INTENT(IN) :: l_ice ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: dtime
  REAL ( KIND= wp ) , INTENT(IN) :: ml_depth
  REAL ( KIND= wp ) , INTENT(IN) :: min_ice_depth
  REAL ( KIND= wp ) :: t_star
  REAL ( KIND= wp ) :: t_freeze
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_lwtr,hflx,hflx_resid,hcap,l_ice_formation,l_ice_melting,l_ice)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( hflx_resid , 1 ) , 1
   IF ( l_ice ( jsb_hori ) ) THEN
!cdir CLAW: RETURN statement transformed for parallel region
   ELSE
    IF ( l_ice_melting ( jsb_hori ) ) THEN
     t_lwtr ( jsb_hori ) = t_lwtr ( jsb_hori ) + dtime * hflx_resid ( jsb_hori&
      ) / hcap ( jsb_hori )
     hflx_resid ( jsb_hori ) = 0._wp
     l_ice_melting ( jsb_hori ) = .FALSE.
    ELSE
     hcap ( jsb_hori ) = rhoh2o * clw * ml_depth
     t_freeze = tmelt - min_ice_depth * rhoi * alf / hcap ( jsb_hori )
     t_star = t_lwtr ( jsb_hori ) + dtime * ( hflx ( jsb_hori ) + hflx_resid (&
      jsb_hori ) ) / hcap ( jsb_hori )
     IF ( t_star >= tmelt ) THEN
      t_lwtr ( jsb_hori ) = t_star
      hflx_resid ( jsb_hori ) = 0._wp
      l_ice_formation ( jsb_hori ) = .FALSE.
     ELSE
      t_lwtr ( jsb_hori ) = tmelt
      hflx_resid ( jsb_hori ) = ( t_star - tmelt ) * hcap ( jsb_hori ) / dtime
      l_ice_formation ( jsb_hori ) = t_star <= t_freeze
     END IF
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_lake_wtr_temp

 SUBROUTINE calc_lice_thickness ( hi , fi , hflx_form_ice , lhflx ,&
  hflx_cond_ice , hflx_melt_lice , l_ice_formation , l_ice_melting , fs ,&
  hi_min , dtime )
  USE mo_jsb_physical_constants , ONLY: rhoi , alf , als

  REAL ( KIND= wp ) , INTENT(INOUT) :: hi ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: fi ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: hflx_form_ice ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: lhflx ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: hflx_cond_ice ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: hflx_melt_lice ( : )
  LOGICAL , INTENT(IN) :: l_ice_formation ( : )
  LOGICAL , INTENT(INOUT) :: l_ice_melting ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: fs ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: hi_min
  REAL ( KIND= wp ) , INTENT(IN) :: dtime
  INTEGER :: jsb_hori

!$acc data &
!$acc present(hi,fi,hflx_form_ice,lhflx,hflx_cond_ice,hflx_melt_lice &
!$acc ,l_ice_formation,l_ice_melting,fs)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( hi , 1 ) , 1
   IF ( l_ice_formation ( jsb_hori ) ) THEN
    hi ( jsb_hori ) = ( - hflx_form_ice ( jsb_hori ) * dtime / ( rhoi * alf )&
     )
    fi ( jsb_hori ) = 1._wp
    hflx_form_ice ( jsb_hori ) = 0._wp
   ELSE
    IF ( hi ( jsb_hori ) >= hi_min ) THEN
     hi ( jsb_hori ) = hi ( jsb_hori ) + dtime / rhoi * ( ( hflx_cond_ice (&
      jsb_hori ) - hflx_melt_lice ( jsb_hori ) ) / alf + ( 1._wp - fs (&
      jsb_hori ) ) * lhflx ( jsb_hori ) / als )
     IF ( hi ( jsb_hori ) >= hi_min ) THEN
      hflx_melt_lice ( jsb_hori ) = 0._wp
      fi ( jsb_hori ) = 1._wp
      l_ice_melting ( jsb_hori ) = .FALSE.
     ELSE
      IF ( hi ( jsb_hori ) <= 0._wp ) THEN
       hflx_melt_lice ( jsb_hori ) = ( - rhoi * alf * hi ( jsb_hori ) / dtime&
        )
       hi ( jsb_hori ) = 0._wp
       fi ( jsb_hori ) = 0._wp
       l_ice_melting ( jsb_hori ) = .TRUE.
      ELSE
       hflx_melt_lice ( jsb_hori ) = rhoi * alf * ( hi_min - hi ( jsb_hori ) )&
        / dtime
       hi ( jsb_hori ) = hi_min
       fi ( jsb_hori ) = 1._wp
       l_ice_melting ( jsb_hori ) = .FALSE.
      END IF
     END IF
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_lice_thickness

 SUBROUTINE calc_lake_ice_temp ( t_lice , hcap_ice , hflx_cond_ice ,&
  hflx_melt_snow , hflx_melt_lice , hflx , lhflx , snowfall , hi , hs , fs ,&
  hi_min , dtime )
  USE mo_jsb_physical_constants , ONLY: tmelt , rhoh2o , rhoi , rhoi , rhos ,&
   als , alf , ki , ks , ci , cs

  REAL ( KIND= wp ) , INTENT(INOUT) :: t_lice ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: hcap_ice ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: hflx_cond_ice ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: hflx_melt_snow ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: hflx_melt_lice ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: hflx ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: lhflx ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: snowfall ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: hi ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: hs ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: fs ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: hi_min
  REAL ( KIND= wp ) , INTENT(IN) :: dtime
  REAL ( KIND= wp ) :: hi_eff
  REAL ( KIND= wp ) :: hflx_melt
  REAL ( KIND= wp ) :: delta_hs
  REAL ( KIND= wp ) :: delta_t_lice
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_lice,hcap_ice,hflx_cond_ice,hflx_melt_snow,hflx_melt_lice &
!$acc ,hflx,lhflx,snowfall,hi,hs,fs)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( hi , 1 ) , 1
   hcap_ice ( jsb_hori ) = rhoi * ci * hi_min + rhoh2o * cs * hs ( jsb_hori )
   hi_eff = hi ( jsb_hori ) + ki / ks * ( rhoh2o / rhos ) * hs ( jsb_hori )
   IF ( hi ( jsb_hori ) < hi_min ) THEN
    t_lice ( jsb_hori ) = tmelt
    hflx_cond_ice ( jsb_hori ) = 0._wp
    hflx_melt_snow ( jsb_hori ) = 0._wp
    hflx_melt_lice ( jsb_hori ) = 0._wp
    hs ( jsb_hori ) = 0._wp
!cdir CLAW: RETURN statement transformed for parallel region
   ELSE
    hs ( jsb_hori ) = hs ( jsb_hori ) + dtime / rhoh2o * ( snowfall ( jsb_hori&
     ) + fs ( jsb_hori ) * lhflx ( jsb_hori ) / als )
    hs ( jsb_hori ) = max ( hs ( jsb_hori ) , 0._wp )
    t_lice ( jsb_hori ) = ( hflx ( jsb_hori ) + hflx_melt_lice ( jsb_hori ) +&
     hcap_ice ( jsb_hori ) * t_lice ( jsb_hori ) / dtime + ki * tmelt / hi_eff&
     ) / ( hcap_ice ( jsb_hori ) / dtime + ki / hi_eff )
    hflx_melt_snow ( jsb_hori ) = 0._wp
    hflx_melt_lice ( jsb_hori ) = 0._wp
    IF ( t_lice ( jsb_hori ) > tmelt .AND. hs ( jsb_hori ) > 0._wp ) THEN
     hflx_melt = ( hcap_ice ( jsb_hori ) / dtime + ki / hi_eff ) * ( t_lice (&
      jsb_hori ) - tmelt )
     delta_hs = min ( hs ( jsb_hori ) , dtime * hflx_melt / alf / rhoh2o )
     hflx_melt_snow ( jsb_hori ) = delta_hs * rhoh2o * alf / dtime
     delta_t_lice = hflx_melt_snow ( jsb_hori ) / ( hcap_ice ( jsb_hori ) /&
      dtime + ki / hi_eff )
     t_lice ( jsb_hori ) = t_lice ( jsb_hori ) - delta_t_lice
     hs ( jsb_hori ) = hs ( jsb_hori ) - delta_hs
     hi_eff = hi ( jsb_hori ) + ki / ks * ( rhoh2o / rhos ) * hs ( jsb_hori )
     hcap_ice ( jsb_hori ) = rhoi * ci * hi_min + rhoh2o * cs * hs ( jsb_hori&
      )
    END IF
    IF ( t_lice ( jsb_hori ) > tmelt ) THEN
     hflx_melt_lice ( jsb_hori ) = ( hcap_ice ( jsb_hori ) / dtime + ki /&
      hi_eff ) * ( t_lice ( jsb_hori ) - tmelt )
     t_lice ( jsb_hori ) = tmelt
    END IF
    hflx_cond_ice ( jsb_hori ) = ( - ki * ( t_lice ( jsb_hori ) - tmelt ) /&
     hi_eff )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_lake_ice_temp

 SUBROUTINE update_surface_fluxes_lake ( tile , options )
  USE mo_phy_schemes , ONLY: heat_transfer_coef
  USE mo_jsb_physical_constants , ONLY: alv , als

  CLASS ( t_jsb_tile_abstract ) , INTENT(INOUT) :: tile
  TYPE ( t_jsb_task_options ) , INTENT(IN) :: options
  TYPE ( t_seb_memory ) , POINTER :: seb__mem
  TYPE ( t_hydro_memory ) , POINTER :: hydro__mem
  TYPE ( t_a2l_memory ) , POINTER :: a2l__mem
  REAL ( KIND= wp ) , POINTER :: sensible_hflx ( : )
  REAL ( KIND= wp ) , POINTER :: latent_hflx ( : )
  REAL ( KIND= wp ) , POINTER :: s_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: t_acoef_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: t_bcoef_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: drag_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: latent_hflx_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: sensible_hflx_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: evapo_wtr ( : )
  REAL ( KIND= wp ) , POINTER :: s_ice ( : )
  REAL ( KIND= wp ) , POINTER :: t_acoef_ice ( : )
  REAL ( KIND= wp ) , POINTER :: t_bcoef_ice ( : )
  REAL ( KIND= wp ) , POINTER :: drag_ice ( : )
  REAL ( KIND= wp ) , POINTER :: latent_hflx_ice ( : )
  REAL ( KIND= wp ) , POINTER :: sensible_hflx_ice ( : )
  REAL ( KIND= wp ) , POINTER :: evapo_ice ( : )
  REAL ( KIND= wp ) , POINTER :: fract_lice ( : )
  REAL ( KIND= wp ) :: s_air ( 1 : options % nc )
  REAL ( KIND= wp ) :: heat_tcoef ( 1 : options % nc )
  INTEGER :: iblk
  INTEGER :: ics
  INTEGER :: ice
  REAL ( KIND= wp ) :: steplen
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_seb_lake" //&
  ":update_surface_fluxes_lake"
  INTEGER :: claw_induction_2

  iblk = options % iblk
  ics = options % ics
  ice = options % ice
  steplen = options % steplen
  IF ( ( .NOT. tile % is_process_active ( seb_ ) ) ) THEN
   RETURN
  END IF
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_lake" //":update_surface_fluxes_lake" ) ,&
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
  t_acoef_wtr => a2l__mem % t_acoef_wtr % ptr ( ics : ice , iblk )
  t_bcoef_wtr => a2l__mem % t_bcoef_wtr % ptr ( ics : ice , iblk )
  t_acoef_ice => a2l__mem % t_acoef_ice % ptr ( ics : ice , iblk )
  t_bcoef_ice => a2l__mem % t_bcoef_ice % ptr ( ics : ice , iblk )
  drag_wtr => a2l__mem % drag_wtr % ptr ( ics : ice , iblk )
  drag_ice => a2l__mem % drag_ice % ptr ( ics : ice , iblk )
  s_wtr => seb__mem % s_wtr % ptr ( ics : ice , iblk )
  s_ice => seb__mem % s_ice % ptr ( ics : ice , iblk )
  sensible_hflx => seb__mem % sensible_hflx % ptr ( ics : ice , iblk )
  sensible_hflx_wtr => seb__mem % sensible_hflx_wtr % ptr ( ics : ice , iblk )
  latent_hflx => seb__mem % latent_hflx % ptr ( ics : ice , iblk )
  latent_hflx_wtr => seb__mem % latent_hflx_wtr % ptr ( ics : ice , iblk )
  sensible_hflx_ice => seb__mem % sensible_hflx_ice % ptr ( ics : ice , iblk )
  latent_hflx_ice => seb__mem % latent_hflx_ice % ptr ( ics : ice , iblk )
  evapo_wtr => hydro__mem % evapo_wtr % ptr ( ics : ice , iblk )
  evapo_ice => hydro__mem % evapo_ice % ptr ( ics : ice , iblk )
  fract_lice => seb__mem % fract_lice % ptr ( ics : ice , iblk )
!$acc data create(s_air, heat_tcoef)
!$acc data &
!$acc present(s_air(:),sensible_hflx(:),drag_wtr(:),t_bcoef_wtr(:) &
!$acc ,t_acoef_ice(:),sensible_hflx_ice(:),evapo_wtr(:),s_wtr(:) &
!$acc ,latent_hflx_ice(:),latent_hflx(:),t_bcoef_ice(:),drag_ice(:) &
!$acc ,t_acoef_wtr(:),evapo_ice(:),heat_tcoef(:),sensible_hflx_wtr(:) &
!$acc ,s_ice(:),latent_hflx_wtr(:),fract_lice(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_2 = 1 , size ( s_air , 1 )
   s_air ( claw_induction_2 ) = t_acoef_wtr ( claw_induction_2 ) * s_wtr (&
    claw_induction_2 ) + t_bcoef_wtr ( claw_induction_2 )
   heat_tcoef ( claw_induction_2 ) = heat_transfer_coef ( drag_wtr (&
    claw_induction_2 ) , steplen )
   sensible_hflx_wtr ( claw_induction_2 ) = heat_tcoef ( claw_induction_2 ) *&
    ( s_air ( claw_induction_2 ) - s_wtr ( claw_induction_2 ) )
   latent_hflx_wtr ( claw_induction_2 ) = alv * evapo_wtr ( claw_induction_2 )
   s_air ( claw_induction_2 ) = t_acoef_ice ( claw_induction_2 ) * s_ice (&
    claw_induction_2 ) + t_bcoef_ice ( claw_induction_2 )
   heat_tcoef ( claw_induction_2 ) = heat_transfer_coef ( drag_ice (&
    claw_induction_2 ) , steplen )
   sensible_hflx_ice ( claw_induction_2 ) = heat_tcoef ( claw_induction_2 ) *&
    ( s_air ( claw_induction_2 ) - s_ice ( claw_induction_2 ) )
   latent_hflx_ice ( claw_induction_2 ) = als * evapo_ice ( claw_induction_2 )
   sensible_hflx ( claw_induction_2 ) = ( 1._wp - fract_lice (&
    claw_induction_2 ) ) * sensible_hflx_wtr ( claw_induction_2 ) + fract_lice&
    ( claw_induction_2 ) * sensible_hflx_ice ( claw_induction_2 )
   latent_hflx ( claw_induction_2 ) = ( 1._wp - fract_lice ( claw_induction_2&
    ) ) * latent_hflx_wtr ( claw_induction_2 ) + fract_lice ( claw_induction_2&
    ) * latent_hflx_ice ( claw_induction_2 )
  END DO
!$acc end parallel
!$acc end data
!$acc end data
  IF ( debug_on ( ) .AND. iblk == 1 ) THEN
   CALL message ( trim ("mo_seb_lake" //":update_surface_fluxes_lake" ) ,&
   "Finished." )
  END IF
 END SUBROUTINE update_surface_fluxes_lake

END MODULE mo_seb_lake

