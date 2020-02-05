MODULE mo_hydro_process
 USE mo_kind , ONLY: wp
 USE mo_exception , ONLY: finish , message_text
 PRIVATE :: soilhyd
 PRIVATE :: digest_evapotrans
 PRIVATE :: digest_infilt
 PRIVATE :: get_snow_fract_noforest
 PRIVATE :: get_canopy_cond_unstressed_simple
 PRIVATE :: get_canopy_cond_stressed_simple
 PUBLIC :: calc_surface_hydrology_land
 PUBLIC :: calc_surface_hydrology_glacier
 PUBLIC :: calc_soil_hydrology
 PUBLIC :: calc_wskin_fractions_lice
 PUBLIC :: calc_wskin_fractions_veg
 PUBLIC :: calc_wskin_fractions_bare
 PUBLIC :: get_canopy_conductance
 PUBLIC :: get_water_stress_factor
 REAL ( KIND= wp ) , PRIVATE , PARAMETER :: zwdmin = 0.05_wp
 REAL ( KIND= wp ) , PRIVATE , PARAMETER :: zwdcri = 0.9_wp
 REAL ( KIND= wp ) , PRIVATE , PARAMETER :: zdrmin = 0.001_wp / ( 3600._wp *&
  1000._wp )
 REAL ( KIND= wp ) , PRIVATE , PARAMETER :: zdrmax = 0.1_wp / ( 3600._wp *&
  1000._wp )
 REAL ( KIND= wp ) , PRIVATE , PARAMETER :: zdrexp = 1.5_wp
 REAL ( KIND= wp ) , PRIVATE , PARAMETER :: zfcmin = 1.e-10_wp
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_hydro_process"
 INTERFACE get_canopy_conductance
  MODULE PROCEDURE get_canopy_cond_unstressed_simple
  MODULE PROCEDURE get_canopy_cond_stressed_simple
 END INTERFACE

CONTAINS
 SUBROUTINE calc_surface_hydrology_land ( lstart , delta_time , l_dynsnow ,&
  t_srf_unfilt , wind_10m , t_air , skinres_canopy_max , skinres_max ,&
  sfract_srf , wfract_srf , hcap_grnd , evapotrans , evapopot , trans , rain ,&
  snow , skinres , snow_soil , canopy_snow , snow_soil_dens , q_snocpymlt ,&
  snow_accum , snow_melt_soil , evapotrans_soil , evapo_skin , evapo_snow ,&
  water_to_soil , evapo_deficit )
  USE mo_jsb_physical_constants , ONLY: tmelt , rhoh2o , alf , dens_snow_min ,&
   dens_snow_max , dens_snow
  USE mo_hydro_constants , ONLY: interceptionefficiency
  USE mo_sse_constants , ONLY: snow_depth_min

  LOGICAL , INTENT(IN) :: lstart
  LOGICAL , INTENT(IN) :: l_dynsnow
  REAL ( KIND= wp ) , INTENT(IN) :: delta_time
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf_unfilt ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: wind_10m ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_air ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: skinres_canopy_max ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: skinres_max ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: sfract_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: wfract_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: hcap_grnd ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: evapotrans ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: evapopot ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: trans ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: rain ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: snow ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: skinres ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: snow_soil ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: canopy_snow ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: snow_soil_dens ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: q_snocpymlt ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: snow_accum ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: snow_melt_soil ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: evapotrans_soil ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: evapo_skin ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: evapo_snow ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: water_to_soil ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: evapo_deficit ( : )
  REAL ( KIND= wp ) :: rain_in_m
  REAL ( KIND= wp ) :: new_snow
  REAL ( KIND= wp ) :: snow_soil_old
  REAL ( KIND= wp ) :: evapotrans_in_m
  REAL ( KIND= wp ) :: evapotrans_soil_in_m
  REAL ( KIND= wp ) :: evapo_soil_in_m
  REAL ( KIND= wp ) :: evapo_skin_in_m
  REAL ( KIND= wp ) :: trans_in_m
  REAL ( KIND= wp ) :: evapo_snow_in_m
  REAL ( KIND= wp ) :: snow_melt_can
  REAL ( KIND= wp ) :: snow_melt_pot_in_m
  REAL ( KIND= wp ) :: new_snow_canopy
  REAL ( KIND= wp ) :: new_snow_soil
  REAL ( KIND= wp ) :: exp_t
  REAL ( KIND= wp ) :: exp_w
  REAL ( KIND= wp ) :: snow_blown
  REAL ( KIND= wp ) :: canopy_pre_snow
  REAL ( KIND= wp ) :: evapotrans_no_snow
  REAL ( KIND= wp ) :: evapo_snow_pot_in_m
  REAL ( KIND= wp ) :: evapo_skin_pot_in_m
  REAL ( KIND= wp ) :: rain_to_skinres
  REAL ( KIND= wp ) :: skinres_pre_rain
  REAL ( KIND= wp ) , PARAMETER :: zc1 = tmelt - 3._wp
  REAL ( KIND= wp ) , PARAMETER :: zc2 = 1.87e5_wp
  REAL ( KIND= wp ) , PARAMETER :: zc3 = 1.56e5_wp
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_srf_unfilt,wind_10m,t_air,skinres_canopy_max,skinres_max &
!$acc ,sfract_srf,wfract_srf,hcap_grnd,evapotrans,evapopot,trans,rain,snow &
!$acc ,skinres,snow_soil,canopy_snow,snow_soil_dens,q_snocpymlt,snow_accum &
!$acc ,snow_melt_soil,evapotrans_soil,evapo_skin,evapo_snow,water_to_soil &
!$acc ,evapo_deficit)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( snow_melt_soil , 1 ) , 1
   snow_accum ( jsb_hori ) = 0._wp
   snow_soil_old = snow_soil ( jsb_hori )
   snow_melt_soil ( jsb_hori ) = 0._wp
   evapotrans_soil ( jsb_hori ) = 0._wp
   evapo_skin ( jsb_hori ) = 0._wp
   evapo_snow ( jsb_hori ) = 0._wp
   water_to_soil ( jsb_hori ) = 0._wp
   evapo_deficit ( jsb_hori ) = 0._wp
   rain_in_m = rain ( jsb_hori ) * delta_time / rhoh2o
   new_snow = snow ( jsb_hori ) * delta_time / rhoh2o
   evapotrans_in_m = evapotrans ( jsb_hori ) * delta_time / rhoh2o
   trans_in_m = trans ( jsb_hori ) * delta_time / rhoh2o
   evapo_snow_pot_in_m = sfract_srf ( jsb_hori ) * evapopot ( jsb_hori ) *&
    delta_time / rhoh2o
   evapotrans_no_snow = evapotrans_in_m - evapo_snow_pot_in_m
   evapo_skin_pot_in_m = ( 1._wp - sfract_srf ( jsb_hori ) ) * wfract_srf (&
    jsb_hori ) * evapopot ( jsb_hori ) * delta_time / rhoh2o
   new_snow_canopy = min ( new_snow * interceptionefficiency , max (&
    skinres_canopy_max ( jsb_hori ) - canopy_snow ( jsb_hori ) , 0._wp ) )
   new_snow_soil = new_snow - new_snow_canopy
   canopy_pre_snow = canopy_snow ( jsb_hori )
   canopy_snow ( jsb_hori ) = min ( max ( 0._wp , canopy_pre_snow +&
    new_snow_canopy + evapo_snow_pot_in_m ) , skinres_canopy_max ( jsb_hori )&
    )
   evapo_snow_in_m = evapo_snow_pot_in_m - ( canopy_snow ( jsb_hori ) -&
    canopy_pre_snow - new_snow_canopy )
   exp_t = max ( 0._wp , t_air ( jsb_hori ) - zc1 ) / zc2 * delta_time
   snow_melt_can = canopy_snow ( jsb_hori ) * ( 1._wp - exp ( ( - exp_t ) ) )
   canopy_snow ( jsb_hori ) = canopy_snow ( jsb_hori ) - snow_melt_can
   exp_w = wind_10m ( jsb_hori ) / zc3 * delta_time
   snow_blown = canopy_snow ( jsb_hori ) * ( 1._wp - exp ( ( - exp_w ) ) )
   canopy_snow ( jsb_hori ) = canopy_snow ( jsb_hori ) - snow_blown
   new_snow_soil = new_snow_soil + snow_blown
   q_snocpymlt ( jsb_hori ) = snow_melt_can * rhoh2o * alf / delta_time
   snow_soil ( jsb_hori ) = snow_soil ( jsb_hori ) + new_snow_soil +&
    evapo_snow_in_m
   IF ( snow_soil ( jsb_hori ) < 0._wp ) THEN
    evapotrans_no_snow = evapotrans_no_snow + snow_soil ( jsb_hori )
    evapo_deficit ( jsb_hori ) = snow_soil ( jsb_hori )
    snow_soil ( jsb_hori ) = 0._wp
   ELSE
    evapo_deficit ( jsb_hori ) = 0._wp
   END IF
   IF ( ( .NOT. lstart ) ) THEN
    IF ( t_srf_unfilt ( jsb_hori ) > tmelt ) THEN
     snow_melt_pot_in_m = hcap_grnd ( jsb_hori ) * ( t_srf_unfilt ( jsb_hori )&
      - tmelt ) / ( alf * rhoh2o )
     IF ( snow_soil ( jsb_hori ) > 0._wp ) THEN
      snow_melt_soil ( jsb_hori ) = max ( min ( snow_melt_pot_in_m , snow_soil&
       ( jsb_hori ) ) , 0._wp )
      snow_soil ( jsb_hori ) = snow_soil ( jsb_hori ) - snow_melt_soil (&
       jsb_hori )
     END IF
    END IF
   END IF
   skinres ( jsb_hori ) = skinres ( jsb_hori ) + snow_melt_can
   water_to_soil ( jsb_hori ) = snow_melt_soil ( jsb_hori ) + max ( 0._wp ,&
    skinres ( jsb_hori ) - skinres_max ( jsb_hori ) )
   skinres ( jsb_hori ) = min ( skinres_max ( jsb_hori ) , skinres ( jsb_hori&
    ) )
   snow_accum ( jsb_hori ) = new_snow + evapo_snow_pot_in_m - snow_melt_soil (&
    jsb_hori ) - snow_melt_can
   rain_to_skinres = min ( rain_in_m * interceptionefficiency , max (&
    skinres_max ( jsb_hori ) - skinres ( jsb_hori ) , 0._wp ) )
   water_to_soil ( jsb_hori ) = water_to_soil ( jsb_hori ) + max ( rain_in_m -&
    rain_to_skinres , 0._wp )
   skinres_pre_rain = skinres ( jsb_hori )
   skinres ( jsb_hori ) = min ( max ( 0._wp , skinres_pre_rain +&
    rain_to_skinres + evapo_skin_pot_in_m ) , skinres_max ( jsb_hori ) )
   evapo_skin_in_m = skinres ( jsb_hori ) - ( skinres_pre_rain +&
    rain_to_skinres )
   evapotrans_soil_in_m = evapotrans_no_snow - evapo_skin_in_m
   IF ( evapo_skin_pot_in_m < 0._wp ) THEN
    evapo_deficit ( jsb_hori ) = evapo_deficit ( jsb_hori ) +&
     evapo_skin_pot_in_m - evapo_skin_in_m
   END IF
   evapo_soil_in_m = evapotrans_soil_in_m - trans_in_m
   IF ( evapo_soil_in_m > 0._wp ) THEN
    water_to_soil ( jsb_hori ) = water_to_soil ( jsb_hori ) + evapo_soil_in_m
   END IF
   IF ( trans_in_m > 0._wp ) THEN
    water_to_soil ( jsb_hori ) = water_to_soil ( jsb_hori ) + trans_in_m
   END IF
   IF ( l_dynsnow ) THEN
    IF ( snow_soil_old * rhoh2o / snow_soil_dens ( jsb_hori ) > snow_depth_min&
     ) THEN
     snow_soil_dens ( jsb_hori ) = dens_snow_max - ( dens_snow_max -&
      snow_soil_dens ( jsb_hori ) ) * exp ( ( - 0.01_wp * delta_time /&
      3600._wp ) )
     IF ( snow_soil ( jsb_hori ) > snow_soil_old ) THEN
      snow_soil_dens ( jsb_hori ) = ( snow_soil_dens ( jsb_hori ) *&
       snow_soil_old + dens_snow_min * ( snow_soil ( jsb_hori ) -&
       snow_soil_old ) ) / snow_soil ( jsb_hori )
     END IF
    ELSE
     snow_soil_dens ( jsb_hori ) = dens_snow_min
    END IF
   ELSE
    snow_soil_dens ( jsb_hori ) = dens_snow
   END IF
   snow_melt_soil ( jsb_hori ) = snow_melt_soil ( jsb_hori ) * rhoh2o /&
    delta_time
   evapotrans_soil ( jsb_hori ) = evapotrans_soil_in_m * rhoh2o / delta_time
   evapo_skin ( jsb_hori ) = evapo_skin_in_m * rhoh2o / delta_time
   evapo_snow ( jsb_hori ) = sfract_srf ( jsb_hori ) * evapopot ( jsb_hori )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_surface_hydrology_land

 SUBROUTINE calc_surface_hydrology_glacier ( lstart , delta_time ,&
  t_srf_unfilt , sfract_srf , hcap_grnd , evapotrans , evapopot , rain , snow&
  , glacier_depth , q_snocpymlt , glacier_melt , glacier_runoff , pme_glacier&
  )
  USE mo_jsb_physical_constants , ONLY: tmelt , rhoh2o , alf

  LOGICAL , INTENT(IN) :: lstart
  REAL ( KIND= wp ) , INTENT(IN) :: delta_time
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf_unfilt ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: sfract_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: hcap_grnd ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: evapotrans ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: evapopot ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: rain ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: snow ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: glacier_depth ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: q_snocpymlt ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: glacier_melt ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: glacier_runoff ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: pme_glacier ( : )
  REAL ( KIND= wp ) :: rain_in_m
  REAL ( KIND= wp ) :: new_snow
  REAL ( KIND= wp ) :: evapotrans_in_m
  REAL ( KIND= wp ) :: snow_melt_pot_in_m
  REAL ( KIND= wp ) :: evapo_snow_pot_in_m
  REAL ( KIND= wp ) :: pme_glacier_in_m
  REAL ( KIND= wp ) :: glacier_runoff_in_m
  REAL ( KIND= wp ) :: glacier_melt_in_m
  REAL ( KIND= wp ) , PARAMETER :: zc1 = tmelt - 3._wp
  REAL ( KIND= wp ) , PARAMETER :: zc2 = 1.87e5_wp
  REAL ( KIND= wp ) , PARAMETER :: zc3 = 1.56e5_wp
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_srf_unfilt,sfract_srf,hcap_grnd,evapotrans,evapopot,rain,snow &
!$acc ,glacier_depth,q_snocpymlt,glacier_melt,glacier_runoff,pme_glacier)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( glacier_depth , 1 ) , 1
   rain_in_m = rain ( jsb_hori ) * delta_time / rhoh2o
   new_snow = snow ( jsb_hori ) * delta_time / rhoh2o
   evapotrans_in_m = evapotrans ( jsb_hori ) * delta_time / rhoh2o
   evapo_snow_pot_in_m = sfract_srf ( jsb_hori ) * evapopot ( jsb_hori ) *&
    delta_time / rhoh2o
   glacier_depth ( jsb_hori ) = glacier_depth ( jsb_hori ) + new_snow +&
    evapo_snow_pot_in_m
   pme_glacier_in_m = rain_in_m + new_snow + evapotrans_in_m
   glacier_runoff_in_m = rain_in_m
   glacier_melt_in_m = 0._wp
   IF ( ( .NOT. lstart ) ) THEN
    IF ( t_srf_unfilt ( jsb_hori ) > tmelt ) THEN
     snow_melt_pot_in_m = hcap_grnd ( jsb_hori ) * ( t_srf_unfilt ( jsb_hori )&
      - tmelt ) / ( alf * rhoh2o )
     glacier_melt_in_m = snow_melt_pot_in_m
     glacier_depth ( jsb_hori ) = glacier_depth ( jsb_hori ) -&
      glacier_melt_in_m
     glacier_runoff_in_m = glacier_runoff_in_m + glacier_melt_in_m
    END IF
   END IF
   q_snocpymlt ( jsb_hori ) = 0._wp
   glacier_melt ( jsb_hori ) = glacier_melt_in_m * rhoh2o / delta_time
   glacier_runoff ( jsb_hori ) = glacier_runoff_in_m * rhoh2o / delta_time
   pme_glacier ( jsb_hori ) = pme_glacier_in_m * rhoh2o / delta_time
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_surface_hydrology_glacier

 SUBROUTINE calc_soil_hydrology ( nidx , nsoil , dz , nlat , delta_time ,&
  ltpe_closed , ltpe_open , soil_depth_l , root_depth_l , fract_org_l ,&
  hyd_cond_sat_l , matrix_pot_l , bclapp_l , pore_size_index_l ,&
  vol_porosity_l , vol_field_cap_l , vol_p_wilt_l , ws_max_mineral ,&
  oro_stddev , t_soil , water_to_soil , evapotrans_soil , trans , ws , ws_l ,&
  soil_ice_l , w_soil_overflow , wsat_l , field_cap_l , p_wilt_l , runoff ,&
  drainage )
  USE mo_jsb_physical_constants , ONLY: tmelt , rhoh2o
  USE mo_hydro_constants , ONLY: oro_var_min , oro_var_max
  USE mo_hydro_util , ONLY: get_water_in_root_zone

  INTEGER , INTENT(IN) :: nidx
  INTEGER , INTENT(IN) :: nsoil
  INTEGER , INTENT(IN) :: nlat
  REAL ( KIND= wp ) , INTENT(IN) :: dz ( 1 : nsoil )
  REAL ( KIND= wp ) , INTENT(IN) :: delta_time
  LOGICAL , INTENT(IN) :: ltpe_closed
  LOGICAL , INTENT(IN) :: ltpe_open
  REAL ( KIND= wp ) , INTENT(IN) :: soil_depth_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: root_depth_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: fract_org_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: hyd_cond_sat_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: matrix_pot_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: bclapp_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: pore_size_index_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: vol_porosity_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: vol_field_cap_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: vol_p_wilt_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: ws_max_mineral ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: oro_stddev ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: water_to_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: evapotrans_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: trans ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: ws ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: w_soil_overflow ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: ws_l ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: soil_ice_l ( : , : )
  REAL ( KIND= wp ) , INTENT(OUT) :: runoff ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: drainage ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: wsat_l ( : , : )
  REAL ( KIND= wp ) , INTENT(OUT) :: field_cap_l ( : , : )
  REAL ( KIND= wp ) , INTENT(OUT) :: p_wilt_l ( : , : )
  REAL ( KIND= wp ) :: sigma_0
  REAL ( KIND= wp ) :: sigma_max
  REAL ( KIND= wp ) :: steepness
  REAL ( KIND= wp ) :: zbws
  REAL ( KIND= wp ) :: zb1
  REAL ( KIND= wp ) :: zbm
  REAL ( KIND= wp ) :: zconw1
  REAL ( KIND= wp ) :: zvol
  REAL ( KIND= wp ) :: evapo_soil_in_m ( 1 : nidx )
  REAL ( KIND= wp ) :: trans_in_m ( 1 : nidx )
  REAL ( KIND= wp ) :: ws_rel ( 1 : nidx )
  REAL ( KIND= wp ) :: infilt ( 1 : nidx )
  REAL ( KIND= wp ) :: ws_max ( 1 : nidx )
  REAL ( KIND= wp ) :: ws_min_drain ( 1 : nidx )
  REAL ( KIND= wp ) :: ws_incl_ice ( 1 : nidx )
  REAL ( KIND= wp ) :: ice_rootzone ( 1 : nidx )
  REAL ( KIND= wp ) :: vfc_times_sd_l ( 1 : nidx , 1 : nsoil )
  INTEGER :: jl
  INTEGER :: i
  INTEGER , PARAMETER :: ilog = 0
  INTEGER :: jllog = 0
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1

!$acc data &
!$acc pcreate(evapo_soil_in_m,trans_in_m,infilt,ws_max,ws_min_drain &
!$acc ,ws_incl_ice,ice_rootzone,vfc_times_sd_l)
  sigma_0 = oro_var_min
  sigma_max = oro_var_max * 64._wp / real ( nlat , wp )
!$acc data &
!$acc present(trans_in_m,trans,wsat_l,vol_porosity_l,soil_depth_l,soil_ice_l &
!$acc ,field_cap_l,vol_field_cap_l,p_wilt_l,vol_p_wilt_l)
!$acc parallel
!$acc loop gang
  DO i = 1 , nidx , 1
   trans_in_m ( i ) = trans ( i ) * delta_time / rhoh2o
!$acc loop seq
   DO jl = 1 , nsoil , 1
    wsat_l ( i , jl ) = max ( vol_porosity_l ( i , jl ) * soil_depth_l ( i ,&
     jl ) - soil_ice_l ( i , jl ) , 0._wp )
    field_cap_l ( i , jl ) = max ( vol_field_cap_l ( i , jl ) * soil_depth_l (&
     i , jl ) - soil_ice_l ( i , jl ) , 0._wp )
    IF ( vol_porosity_l ( i , jl ) > 0._wp ) THEN
     p_wilt_l ( i , jl ) = wsat_l ( i , jl ) * vol_p_wilt_l ( i , jl ) /&
      vol_porosity_l ( i , jl )
    ELSE
     p_wilt_l ( i , jl ) = 0._wp
    END IF
   END DO
  END DO
!$acc end parallel
!$acc end data
!$acc data present(vol_field_cap_l(:,:),soil_depth_l(:,:),vfc_times_sd_l(:,:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_0 = 1 , size ( vfc_times_sd_l , 1 )
   DO claw_induction_1 = 1 , size ( vfc_times_sd_l , 2 )
    vfc_times_sd_l ( claw_induction_0 , claw_induction_1 ) = vol_field_cap_l (&
     claw_induction_0 , claw_induction_1 ) * soil_depth_l ( claw_induction_0 ,&
     claw_induction_1 )
   END DO
  END DO
!$acc end parallel
!$acc end data
  CALL get_water_in_root_zone ( vfc_times_sd_l ( : , : ) , soil_depth_l ( : ,&
   : ) , root_depth_l ( : , : ) , ws_max ( : ) )
!$acc data present(fract_org_l,ws_max,ws_max_mineral)
!$acc parallel loop gang vector
  DO i = 1 , nidx , 1
   IF ( sum ( fract_org_l ( i , : ) ) <= epsilon ( 1._wp ) ) THEN
    ws_max ( i ) = ws_max_mineral ( i )
   END IF
  END DO
!$acc end parallel
!$acc end data
  CALL get_water_in_root_zone ( soil_ice_l ( : , : ) , soil_depth_l ( : , : )&
   , root_depth_l ( : , : ) , ice_rootzone ( : ) )
  CALL get_water_in_root_zone ( ws_l ( : , : ) , soil_depth_l ( : , : ) ,&
   root_depth_l ( : , : ) , ws ( : ) )
!$acc data present(evapo_soil_in_m,evapotrans_soil,trans_in_m) &
!$acc present(ws_min_drain,ws_max,oro_stddev,runoff,infilt) &
!$acc present(water_to_soil,ws_incl_ice,ice_rootzone,drainage,t_soil)
!$acc parallel
!$acc loop gang
  DO i = 1 , nidx , 1
   evapo_soil_in_m ( i ) = evapotrans_soil ( i ) * delta_time / rhoh2o -&
    trans_in_m ( i )
  END DO
!$acc end parallel
!$acc parallel
!$acc loop gang vector
  DO i = 1 , nidx , 1
   ws_min_drain ( i ) = zwdmin * ws_max ( i )
   steepness = max ( 0._wp , oro_stddev ( i ) - sigma_0 ) / ( oro_stddev ( i )&
    + sigma_max )
   zbws = max ( min ( steepness , 0.5_wp ) , 0.01_wp )
   zb1 = 1._wp + zbws
   zbm = 1._wp / zb1
   zconw1 = ws_max ( i ) * zb1
   IF ( ltpe_closed ) THEN
    runoff ( i ) = 0._wp
    infilt ( i ) = water_to_soil ( i )
   ELSE
    IF ( t_soil ( i ) < tmelt ) THEN
     runoff ( i ) = water_to_soil ( i )
     infilt ( i ) = 0._wp
    ELSE
     ws_incl_ice ( i ) = ws ( i ) + ice_rootzone ( i )
     IF ( water_to_soil ( i ) > 0._wp .AND. ws_incl_ice ( i ) > ws_min_drain (&
      i ) ) THEN
      ws_rel ( i ) = min ( 1._wp , ws_incl_ice ( i ) / ws_max ( i ) )
      zvol = ( 1._wp - ws_rel ( i ) ) ** ( zbm ) - water_to_soil ( i ) /&
       zconw1
      runoff ( i ) = water_to_soil ( i ) - ( ws_max ( i ) - ws_incl_ice ( i )&
       )
      IF ( zvol > 0._wp ) THEN
       runoff ( i ) = runoff ( i ) + ws_max ( i ) * zvol ** ( zb1 )
      END IF
      runoff ( i ) = max ( min ( runoff ( i ) , water_to_soil ( i ) ) , 0._wp&
       )
      infilt ( i ) = water_to_soil ( i ) - runoff ( i )
     ELSE
      runoff ( i ) = 0._wp
      infilt ( i ) = water_to_soil ( i )
     END IF
    END IF
   END IF
  END DO
!$acc end parallel
!$acc parallel
!$acc loop gang
  DO i = 1 , nidx , 1
   drainage ( i ) = 0._wp
  END DO
!$acc end parallel
!$acc end data
  CALL digest_evapotrans ( nidx , nsoil , root_depth_l , soil_depth_l ,&
   field_cap_l , trans_in_m , evapo_soil_in_m , ws_l , soil_ice_l )
!$acc data present(infilt)
!$acc parallel loop gang vector
  DO i = 1 , nidx , 1
   infilt ( i ) = infilt ( i ) / 2._wp
  END DO
!$acc end data
  CALL digest_infilt ( nidx , nsoil , ltpe_closed , field_cap_l , infilt ,&
   ws_l , drainage , w_soil_overflow )
  CALL soilhyd ( nidx , nsoil , delta_time , 0 , jllog , ws_l , soil_depth_l ,&
   wsat_l , field_cap_l , hyd_cond_sat_l , vol_porosity_l , bclapp_l ,&
   matrix_pot_l , drainage , pore_size_index_l , p_wilt_l , dz , soil_ice_l ,&
   w_soil_overflow , ltpe_closed , ltpe_open )
  CALL digest_infilt ( nidx , nsoil , ltpe_closed , field_cap_l , infilt ,&
   ws_l , drainage , w_soil_overflow )
  CALL get_water_in_root_zone ( ws_l ( : , : ) , soil_depth_l ( : , : ) ,&
   root_depth_l ( : , : ) , ws ( : ) )
!$acc data present(runoff,drainage)
!$acc parallel loop gang vector
  DO i = 1 , nidx , 1
   runoff ( i ) = runoff ( i ) * rhoh2o / delta_time
   drainage ( i ) = drainage ( i ) * rhoh2o / delta_time
  END DO
!$acc end data
!$acc end data
 END SUBROUTINE calc_soil_hydrology

 SUBROUTINE digest_evapotrans ( nidx , nsoil , droot , dsoil , field_cap_l ,&
  trans_in_m , evapo_soil , ws_l , soil_ice )

  INTEGER , INTENT(IN) :: nidx
  INTEGER , INTENT(IN) :: nsoil
  REAL ( KIND= wp ) , INTENT(IN) :: droot ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: dsoil ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: field_cap_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: trans_in_m ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: evapo_soil ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: ws_l ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: soil_ice ( : , : )
  INTEGER :: jk
  INTEGER :: i
  REAL ( KIND= wp ) :: deficit_evapo ( 1 : nidx )
  REAL ( KIND= wp ) :: deficit_trans ( 1 : nidx )
  REAL ( KIND= wp ) :: deficit_l ( 1 : nidx )
  REAL ( KIND= wp ) :: remaining ( 1 : nidx )
  REAL ( KIND= wp ) :: rootfract ( 1 : nidx )
  REAL ( KIND= wp ) :: fixed ( 1 : nidx )
  REAL ( KIND= wp ) :: trans_layer ( 1 : nidx )
  REAL ( KIND= wp ) :: ws ( 1 : nidx )
  REAL ( KIND= wp ) :: root_depth ( 1 : nidx )
  REAL ( KIND= wp ) :: reduce_evapo ( 1 : nidx )
  INTEGER :: max2d ( 1 : 2 )
  REAL ( KIND= wp ) , PARAMETER :: zwilt = 0.35_wp

!$acc data present(droot, dsoil, field_cap_l, trans_in_m, evapo_soil, &
!$acc ws_l, soil_ice) pcreate(deficit_evapo, deficit_trans, deficit_l, &
!$acc remaining, rootfract, fixed, trans_layer, ws, root_depth, reduce_evapo)
!$acc parallel
!$acc loop gang vector
  DO i = 1 , nidx , 1
   deficit_evapo ( i ) = 0._wp
   IF ( evapo_soil ( i ) < 0._wp ) THEN
    ws_l ( i , 1 ) = ws_l ( i , 1 ) + evapo_soil ( i )
    IF ( ws_l ( i , 1 ) < 0._wp ) THEN
     deficit_evapo ( i ) = ws_l ( i , 1 )
     ws_l ( i , 1 ) = 0._wp
    END IF
   END IF
   root_depth ( i ) = sum ( droot ( i , : ) )
   IF ( root_depth ( i ) <= 0._wp .AND. trans_in_m ( i ) < 0._wp ) THEN
    deficit_trans ( i ) = trans_in_m ( i )
   ELSE
    deficit_trans ( i ) = 0._wp
   END IF
!$acc loop seq
   DO jk = 1 , nsoil , 1
    deficit_l ( i ) = 0._wp
    rootfract ( i ) = 0._wp
    fixed ( i ) = 0._wp
    trans_layer ( i ) = 0._wp
    remaining ( i ) = 0._wp
    IF ( root_depth ( i ) > 0._wp .AND. trans_in_m ( i ) < 0._wp ) THEN
     trans_layer ( i ) = trans_in_m ( i ) * droot ( i , jk ) / root_depth ( i&
      )
     IF ( droot ( i , jk ) > 0._wp ) THEN
      rootfract ( i ) = droot ( i , jk ) / dsoil ( i , jk )
      fixed ( i ) = field_cap_l ( i , jk ) * zwilt * rootfract ( i )
      IF ( ws_l ( i , jk ) * rootfract ( i ) >= fixed ( i ) ) THEN
       remaining ( i ) = ws_l ( i , jk ) * rootfract ( i ) + trans_layer ( i )
       IF ( remaining ( i ) < fixed ( i ) ) THEN
        deficit_l ( i ) = remaining ( i ) - fixed ( i )
        ws_l ( i , jk ) = fixed ( i ) + ws_l ( i , jk ) * ( 1._wp - rootfract&
         ( i ) )
       ELSE
        deficit_l ( i ) = 0._wp
        ws_l ( i , jk ) = remaining ( i ) + ws_l ( i , jk ) * ( 1._wp -&
         rootfract ( i ) )
       END IF
      ELSE
       deficit_l ( i ) = trans_layer ( i )
      END IF
      deficit_trans ( i ) = deficit_trans ( i ) + deficit_l ( i )
     END IF
    END IF
   END DO
   reduce_evapo ( i ) = deficit_evapo ( i ) + deficit_trans ( i )
   ws ( i ) = sum ( ws_l ( i , : ) ) + sum ( soil_ice ( i , : ) )
!$acc loop seq
   DO jk = 1 , nsoil , 1
    IF ( reduce_evapo ( i ) < 0._wp .AND. ws ( i ) > 0._wp ) THEN
     ws_l ( i , jk ) = max ( 0._wp , ws_l ( i , jk ) + reduce_evapo ( i ) *&
      ws_l ( i , jk ) / ws ( i ) )
     soil_ice ( i , jk ) = max ( 0._wp , soil_ice ( i , jk ) + reduce_evapo (&
      i ) * soil_ice ( i , jk ) / ws ( i ) )
    END IF
   END DO
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE digest_evapotrans

 SUBROUTINE digest_infilt ( nidx , nsoil , ltpe_closed , field_cap_l , infilt&
  , ws_l , drainage , w_soil_overflow )

  INTEGER , INTENT(IN) :: nidx
  INTEGER , INTENT(IN) :: nsoil
  LOGICAL , INTENT(IN) :: ltpe_closed
  REAL ( KIND= wp ) , INTENT(IN) :: field_cap_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: infilt ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: ws_l ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: drainage ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: w_soil_overflow ( : )
  INTEGER :: jk
  INTEGER :: i
  REAL ( KIND= wp ) :: excess ( 1 : nidx )

!$acc data &
!$acc present(field_cap_l, infilt, ws_l, drainage, w_soil_overflow) &
!$acc create(excess)
!$acc parallel
!$acc loop gang vector
  DO i = 1 , nidx , 1
   excess ( i ) = infilt ( i )
!$acc loop seq
   DO jk = 1 , nsoil , 1
    ws_l ( i , jk ) = ws_l ( i , jk ) + excess ( i )
    IF ( ( .NOT. ltpe_closed ) ) THEN
     IF ( ws_l ( i , jk ) > field_cap_l ( i , jk ) ) THEN
      excess ( i ) = ws_l ( i , jk ) - field_cap_l ( i , jk )
      ws_l ( i , jk ) = field_cap_l ( i , jk )
     ELSE
      excess ( i ) = 0._wp
     END IF
    ELSE
     IF ( ws_l ( i , jk ) > field_cap_l ( i , jk ) ) THEN
      w_soil_overflow ( i ) = w_soil_overflow ( i ) + ws_l ( i , jk ) -&
       field_cap_l ( i , jk )
      ws_l ( i , jk ) = field_cap_l ( i , jk )
     END IF
     excess ( i ) = 0._wp
    END IF
   END DO
   IF ( excess ( i ) > 0._wp ) THEN
    drainage ( i ) = drainage ( i ) + excess ( i )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE digest_infilt

 SUBROUTINE soilhyd ( nidx , nsoil , delta_time , ilog , jllog , ws_l , dsoil&
  , wsat_l , field_cap_l , hyd_cond_sat_l , vol_porosity_l , bclapp_l ,&
  matrix_pot_l , drainage , pore_size_index_l , p_wilt_l , dz , soil_ice ,&
  w_soil_overflow , ltpe_closed , ltpe_open )

  INTEGER , INTENT(IN) :: nidx
  INTEGER , INTENT(IN) :: nsoil
  REAL ( KIND= wp ) , INTENT(IN) :: delta_time
  INTEGER , INTENT(IN) :: ilog
  INTEGER , INTENT(IN) :: jllog
  REAL ( KIND= wp ) , INTENT(IN) :: dsoil ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: wsat_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: field_cap_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: hyd_cond_sat_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: vol_porosity_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: bclapp_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: matrix_pot_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: pore_size_index_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: p_wilt_l ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: dz ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: ws_l ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: drainage ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: w_soil_overflow ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: soil_ice ( : , : )
  LOGICAL , INTENT(IN) :: ltpe_closed
  LOGICAL , INTENT(IN) :: ltpe_open
  INTEGER , PARAMETER :: iredu = 0
  REAL ( KIND= wp ) , PARAMETER :: credu = 1.0_wp
  INTEGER :: i
  INTEGER :: j
  REAL ( KIND= wp ) :: percolation ( 1 : nidx , 0 : nsoil )
  REAL ( KIND= wp ) :: diffus_l ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: ws_min_drain ( 1 : nidx )
  REAL ( KIND= wp ) :: ws_rel ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: ws_crit ( 1 : nidx )
  REAL ( KIND= wp ) :: drain_crit ( 1 : nidx )
  REAL ( KIND= wp ) :: ws_tmp ( 1 : nidx )
  REAL ( KIND= wp ) :: drain_l ( 1 : nidx )
  LOGICAL :: soil_above ( 1 : nidx )
  LOGICAL :: soil_below ( 1 : nidx )
  REAL ( KIND= wp ) :: zda ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: zdb ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: zdc ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: ztri ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: diffusivity ( 1 : nidx , 1 : nsoil - 1 )
  REAL ( KIND= wp ) :: zkredu ( 1 : nsoil )
  REAL ( KIND= wp ) :: zdiflog ( 1 : 10 )
  REAL ( KIND= wp ) :: wslog ( 1 : nsoil )
  REAL ( KIND= wp ) :: remoist ( 1 : nidx )
  REAL ( KIND= wp ) :: zmvg_l ( 1 : nidx , 1 : nsoil )

!$acc data &
!$acc present(dsoil, wsat_l, field_cap_l, hyd_cond_sat_l, vol_porosity_l, bclapp_l, matrix_pot_l, &
!$acc pore_size_index_l, p_wilt_l, dz, ws_l, drainage, w_soil_overflow, soil_ice) &
!$acc pcreate(diffus_l, ws_min_drain, ws_rel, ws_crit, drain_crit, ws_tmp, drain_l, soil_above, &
!$acc soil_below, zda, zdb, zdc, ztri, diffusivity, zkredu, zdiflog, remoist, zmvg_l, &
!$acc percolation, wslog)
!$acc parallel loop gang vector collapse(2)
  DO j = 1 , nidx , 1
   DO i = 1 , nsoil , 1
    ws_rel ( j , i ) = 0._wp
   END DO
  END DO
  IF ( 0 == 1 ) THEN
!$acc parallel
   zkredu ( 1 ) = ( 1._wp + exp ( ( - credu * dz ( 1 ) ) ) ) / 2._wp
!$acc loop gang vector
   DO i = 2 , nsoil , 1
    zkredu ( i ) = ( exp ( ( - credu * sum ( dz ( 1 : i - 1 ) ) ) ) + exp ( (&
     - credu * sum ( dz ( 1 : i ) ) ) ) ) / 2._wp
   END DO
!$acc end parallel
  ELSE
!$acc parallel loop gang vector
   DO i = 1 , nsoil , 1
    zkredu ( i ) = 1._wp
   END DO
  END IF
!$acc parallel loop gang vector collapse(2)
  DO j = 1 , nidx , 1
   DO i = 1 , nsoil , 1
    zmvg_l ( j , i ) = pore_size_index_l ( j , i ) / ( pore_size_index_l ( j ,&
     i ) + 1._wp )
   END DO
  END DO
!$acc parallel default(none)
!$acc loop gang vector
  DO j = 1 , nidx , 1
   percolation ( j , nsoil ) = 0._wp
!$acc loop seq
   DO i = 1 , nsoil - 1 , 1
    IF ( field_cap_l ( j , i ) > zfcmin .AND. field_cap_l ( j , i ) > p_wilt_l&
     ( j , i ) ) THEN
     ws_rel ( j , i ) = max ( ws_l ( j , i ) - p_wilt_l ( j , i ) , 0._wp ) /&
      ( field_cap_l ( j , i ) - p_wilt_l ( j , i ) )
     ws_rel ( j , i ) = min ( ws_rel ( j , i ) , 1._wp )
     percolation ( j , i ) = hyd_cond_sat_l ( j , i ) * zkredu ( i ) * sqrt (&
      ws_rel ( j , i ) ) * ( 1._wp - ( 1._wp - ws_rel ( j , i ) ** ( 1._wp /&
      zmvg_l ( j , i ) ) ) ** ( zmvg_l ( j , i ) ) ) ** ( 2._wp ) * delta_time
    ELSE
     percolation ( j , i ) = 0._wp
    END IF
   END DO
!$acc loop seq
   DO i = 1 , nsoil , 1
    soil_above ( j ) = i > 1
    IF ( i == nsoil ) THEN
     soil_below ( j ) = .FALSE.
    ELSE
     soil_below ( j ) = dsoil ( j , i + 1 ) > 0._wp
    END IF
    ws_crit ( j ) = field_cap_l ( j , i ) * zwdcri
    IF ( dsoil ( j , i ) > 0._wp .AND. ( .NOT. soil_below ( j ) ) ) THEN
     IF ( ws_l ( j , i ) <= p_wilt_l ( j , i ) ) THEN
      drain_l ( j ) = 0._wp
     ELSE
      IF ( field_cap_l ( j , i ) > p_wilt_l ( j , i ) ) THEN
       drain_l ( j ) = zdrmin * ( ws_l ( j , i ) - p_wilt_l ( j , i ) ) / (&
        field_cap_l ( j , i ) - p_wilt_l ( j , i ) )
      ELSE
       drain_l ( j ) = 0._wp
      END IF
      IF ( ws_l ( j , i ) > ws_crit ( j ) ) THEN
       drain_crit ( j ) = ( zdrmax - zdrmin ) * ( ( ws_l ( j , i ) - ws_crit (&
        j ) ) / ( wsat_l ( j , i ) - ws_crit ( j ) ) ) ** ( zdrexp )
       drain_l ( j ) = drain_l ( j ) + min ( drain_crit ( j ) , ( ws_l ( j , i&
        ) - ws_crit ( j ) ) / delta_time )
      END IF
      drain_l ( j ) = min ( drain_l ( j ) , ( ws_l ( j , i ) - p_wilt_l ( j ,&
       i ) ) / delta_time )
     END IF
     percolation ( j , i ) = drain_l ( j ) * delta_time
    ELSE
     IF ( dsoil ( j , i ) > 0._wp ) THEN
      IF ( ( .NOT. soil_above ( j ) ) ) THEN
       ws_tmp ( j ) = ws_l ( j , i ) - percolation ( j , i )
      ELSE
       ws_tmp ( j ) = ws_l ( j , i ) - percolation ( j , i ) + percolation ( j&
        , i - 1 )
      END IF
      IF ( ws_tmp ( j ) > p_wilt_l ( j , i ) ) THEN
       IF ( field_cap_l ( j , i ) > p_wilt_l ( j , i ) ) THEN
        drain_l ( j ) = zdrmin * ( ws_tmp ( j ) - p_wilt_l ( j , i ) ) / (&
         field_cap_l ( j , i ) - p_wilt_l ( j , i ) )
       ELSE
        drain_l ( j ) = 0._wp
       END IF
       IF ( ws_tmp ( j ) > ws_crit ( j ) ) THEN
        drain_crit ( j ) = ( zdrmax - zdrmin ) * ( ( ws_tmp ( j ) - ws_crit (&
         j ) ) / ( wsat_l ( j , i ) - ws_crit ( j ) ) ) ** ( zdrexp )
        drain_l ( j ) = drain_l ( j ) + min ( drain_crit ( j ) , ( ws_tmp ( j&
         ) - ws_crit ( j ) ) / delta_time )
       END IF
       drain_l ( j ) = min ( drain_l ( j ) * delta_time , ws_tmp ( j ) -&
        p_wilt_l ( j , i ) )
      END IF
     END IF
    END IF
    IF ( ( .NOT. ltpe_closed ) ) THEN
     IF ( dsoil ( j , i ) > 0._wp .AND. soil_below ( j ) .AND. ws_tmp ( j ) >&
      p_wilt_l ( j , i ) ) THEN
      ws_l ( j , i ) = ws_l ( j , i ) - drain_l ( j )
      drainage ( j ) = drainage ( j ) + drain_l ( j )
     END IF
    END IF
   END DO
!$acc loop seq
   DO i = 1 , nsoil , 1
    ws_rel ( j , i ) = 0._wp
    IF ( ws_l ( j , i ) > 0._wp ) THEN
     ws_rel ( j , i ) = ws_l ( j , i ) / dsoil ( j , i ) / ( vol_porosity_l (&
      j , i ) - soil_ice ( j , i ) / dsoil ( j , i ) )
    END IF
    diffus_l ( j , i ) = 0._wp
    IF ( ws_rel ( j , i ) > zfcmin ) THEN
     diffus_l ( j , i ) = bclapp_l ( j , i ) * hyd_cond_sat_l ( j , i ) *&
      matrix_pot_l ( j , i ) / ( ws_l ( j , i ) / dsoil ( j , i ) ) * ws_rel (&
      j , i ) ** ( bclapp_l ( j , i ) + 3._wp )
    END IF
   END DO
!$acc loop seq
   DO i = 1 , nsoil - 1 , 1
    IF ( diffus_l ( j , i + 1 ) == 0._wp .OR. diffus_l ( j , i ) == 0._wp )&
     THEN
     diffusivity ( j , i ) = 0._wp
    ELSE
     diffusivity ( j , i ) = diffus_l ( j , i + 1 ) * diffus_l ( j , i ) / ( (&
      diffus_l ( j , i + 1 ) * dsoil ( j , i ) + diffus_l ( j , i ) * dsoil (&
      j , i + 1 ) ) / ( dsoil ( j , i + 1 ) + dsoil ( j , i ) ) )
    END IF
   END DO
   zda ( j , nsoil ) = 0._wp
   IF ( dsoil ( j , nsoil ) > 0._wp ) THEN
    zdb ( j , nsoil ) = diffusivity ( j , nsoil - 1 ) * delta_time / dsoil ( j&
     , nsoil ) / ( dsoil ( j , nsoil ) + dsoil ( j , nsoil - 1 ) ) * 2._wp
   ELSE
    zdb ( j , nsoil ) = 0._wp
   END IF
!$acc loop seq
   DO i = nsoil - 1 , 2 , (-1)
    IF ( dsoil ( j , i ) > 0._wp ) THEN
     IF ( dsoil ( j , i + 1 ) > 0._wp ) THEN
      zda ( j , i ) = diffusivity ( j , i ) * delta_time / dsoil ( j , i ) / (&
       dsoil ( j , i ) + dsoil ( j , i + 1 ) ) * 2._wp
     ELSE
      zda ( j , i ) = 0._wp
     END IF
     zdb ( j , i ) = diffusivity ( j , i - 1 ) * delta_time / dsoil ( j , i )&
      / ( dsoil ( j , i ) + dsoil ( j , i - 1 ) ) * 2._wp
    ELSE
     zda ( j , i ) = 0._wp
     zdb ( j , i ) = 0._wp
    END IF
   END DO
   zdb ( j , 1 ) = 0._wp
   IF ( dsoil ( j , 1 ) > 0 ) THEN
    IF ( dsoil ( j , 2 ) > 0 ) THEN
     zda ( j , 1 ) = diffusivity ( j , 1 ) * delta_time / dsoil ( j , 1 ) / (&
      dsoil ( j , 1 ) + dsoil ( j , 2 ) ) * 2._wp
    ELSE
     zda ( j , 1 ) = 0._wp
    END IF
   ELSE
    zda ( j , 1 ) = 0._wp
   END IF
   IF ( ilog == 1 .OR. ilog == 3 ) THEN
    DO i = 1 , nsoil , 1
     wslog ( i ) = ws_l ( jllog , i )
    END DO
   END IF
   IF ( dsoil ( j , 1 ) > 0._wp ) THEN
    zdc ( j , 1 ) = zda ( j , 1 ) + zdb ( j , 1 ) + 1._wp
    ws_l ( j , 1 ) = ws_l ( j , 1 ) / zdc ( j , 1 )
   END IF
!$acc loop seq
   DO i = 2 , nsoil , 1
    IF ( dsoil ( j , i ) > 0._wp ) THEN
     ztri ( j , i ) = ( - zda ( j , i - 1 ) / zdc ( j , i - 1 ) )
     zdc ( j , i ) = zda ( j , i ) + zdb ( j , i ) + 1._wp + zdb ( j , i ) *&
      ztri ( j , i )
     ws_l ( j , i ) = ( ws_l ( j , i ) / dsoil ( j , i ) + zdb ( j , i ) *&
      ws_l ( j , i - 1 ) / dsoil ( j , i - 1 ) ) / zdc ( j , i ) * dsoil ( j ,&
      i )
    END IF
   END DO
!$acc loop seq
   DO i = nsoil - 1 , 1 , (-1)
    IF ( dsoil ( j , i + 1 ) > 0._wp ) THEN
     ws_l ( j , i ) = ws_l ( j , i ) - ztri ( j , i + 1 ) * ws_l ( j , i + 1 )&
      * ( dsoil ( j , i ) / dsoil ( j , i + 1 ) )
    END IF
   END DO
   IF ( ilog == 1 .OR. ilog == 3 ) THEN
    DO i = 1 , nsoil , 1
     zdiflog ( i ) = wslog ( i ) - ws_l ( jllog , i )
    END DO
   END IF
   IF ( ( .NOT. ltpe_closed ) ) THEN
!$acc loop seq
    DO i = 1 , nsoil , 1
     IF ( dsoil ( j , i ) > 0._wp .AND. ws_l ( j , i ) > field_cap_l ( j , i )&
      ) THEN
      drainage ( j ) = drainage ( j ) + ws_l ( j , i ) - field_cap_l ( j , i )
      ws_l ( j , i ) = field_cap_l ( j , i )
     END IF
    END DO
   ELSE
!$acc loop seq
    DO i = 1 , nsoil , 1
     IF ( dsoil ( j , i ) > 0._wp ) THEN
      IF ( ws_l ( j , i ) > field_cap_l ( j , i ) ) THEN
       w_soil_overflow ( j ) = w_soil_overflow ( j ) + ws_l ( j , i ) -&
        field_cap_l ( j , i )
       ws_l ( j , i ) = field_cap_l ( j , i )
      ELSE
       remoist ( j ) = min ( w_soil_overflow ( j ) , field_cap_l ( j , i ) -&
        ws_l ( j , i ) )
       ws_l ( j , i ) = ws_l ( j , i ) + remoist ( j )
       w_soil_overflow ( j ) = w_soil_overflow ( j ) - remoist ( j )
      END IF
     END IF
    END DO
   END IF
   IF ( ltpe_open ) THEN
!$acc loop seq
    DO i = 1 , nsoil , 1
     IF ( i == nsoil - 1 .OR. i == nsoil ) THEN
      IF ( dsoil ( j , i ) > 0._wp ) THEN
       ws_l ( j , i ) = max ( ws_l ( j , i ) , field_cap_l ( j , i ) * 0.9_wp&
        )
      END IF
     END IF
    END DO
   END IF
!$acc loop seq
   DO i = 1 , nsoil , 1
    ws_min_drain ( j ) = 0._wp
    IF ( dsoil ( j , i ) > 0._wp .AND. percolation ( j , i ) > 0._wp ) THEN
     ws_min_drain ( j ) = zwdmin * field_cap_l ( j , i )
     IF ( ws_l ( j , i ) - percolation ( j , i ) < ws_min_drain ( j ) ) THEN
      IF ( ws_l ( j , i ) > ws_min_drain ( j ) ) THEN
       percolation ( j , i ) = ws_l ( j , i ) - ws_min_drain ( j )
       ws_l ( j , i ) = ws_min_drain ( j )
      ELSE
       percolation ( j , i ) = 0._wp
      END IF
     ELSE
      ws_l ( j , i ) = ws_l ( j , i ) - percolation ( j , i )
     END IF
    END IF
   END DO
!$acc loop seq
   DO i = 1 , nsoil - 1 , 1
    IF ( dsoil ( j , i + 1 ) > 0._wp ) THEN
     ws_l ( j , i + 1 ) = ws_l ( j , i + 1 ) + percolation ( j , i )
    END IF
    IF ( ( .NOT. ltpe_closed ) ) THEN
     IF ( dsoil ( j , i + 1 ) <= 0._wp ) THEN
      drainage ( j ) = drainage ( j ) + percolation ( j , i )
      percolation ( j , i ) = 0._wp
     END IF
    ELSE
     IF ( dsoil ( j , i + 1 ) <= 0._wp ) THEN
      w_soil_overflow ( j ) = w_soil_overflow ( j ) + percolation ( j , i )
      percolation ( j , i ) = 0._wp
     END IF
    END IF
   END DO
!$acc loop seq
   DO i = nsoil , 2 , (-1)
    IF ( dsoil ( j , i ) > 0._wp .AND. ws_l ( j , i ) > field_cap_l ( j , i )&
     ) THEN
     percolation ( j , i - 1 ) = max ( percolation ( j , i - 1 ) - ws_l ( j ,&
      i ) + field_cap_l ( j , i ) , 0._wp )
     ws_l ( j , i - 1 ) = ws_l ( j , i - 1 ) + ws_l ( j , i ) - field_cap_l (&
      j , i )
     ws_l ( j , i ) = field_cap_l ( j , i )
    END IF
   END DO
   IF ( ( .NOT. ltpe_closed ) ) THEN
    drainage ( j ) = drainage ( j ) + percolation ( j , nsoil )
   ELSE
    w_soil_overflow ( j ) = w_soil_overflow ( j ) + percolation ( j , nsoil )
   END IF
   IF ( ( .NOT. ltpe_closed ) ) THEN
    IF ( ws_l ( j , 1 ) > field_cap_l ( j , 1 ) ) THEN
     drainage ( j ) = drainage ( j ) + ws_l ( j , 1 ) - field_cap_l ( j , 1 )
     ws_l ( j , 1 ) = field_cap_l ( j , 1 )
    END IF
   ELSE
    IF ( ws_l ( j , 1 ) > field_cap_l ( j , 1 ) ) THEN
     w_soil_overflow ( j ) = w_soil_overflow ( j ) + ws_l ( j , 1 ) -&
      field_cap_l ( j , 1 )
     ws_l ( j , 1 ) = field_cap_l ( j , 1 )
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE soilhyd

 SUBROUTINE calc_wskin_fractions_lice ( w_snow_lice , fract_snow_lice )

  REAL ( KIND= wp ) , INTENT(IN) :: w_snow_lice ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: fract_snow_lice ( : )
  INTEGER :: jsb_hori

!$acc data present(w_snow_lice,fract_snow_lice)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( w_snow_lice , 1 ) , 1
   fract_snow_lice ( jsb_hori ) = tanh ( w_snow_lice ( jsb_hori ) * 100._wp )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_wskin_fractions_lice

 SUBROUTINE calc_wskin_fractions_veg ( dtime , steplen , wsr_max , oro_stddev&
  , t_srf_old , press_srf , drag_srf , q_air , w_skin , w_snow_soil ,&
  w_snow_can , fract_snow_can , fract_water , fract_snow_soil )
  USE mo_phy_schemes , ONLY: heat_transfer_coef , qsat
  USE mo_jsb_physical_constants , ONLY: rhoh2o

  REAL ( KIND= wp ) , INTENT(IN) :: dtime
  REAL ( KIND= wp ) , INTENT(IN) :: steplen
  REAL ( KIND= wp ) , INTENT(IN) :: oro_stddev ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf_old ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: press_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: drag_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: q_air ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: wsr_max ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: w_skin ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: w_snow_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: w_snow_can ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: fract_snow_can ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: fract_water ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: fract_snow_soil ( : )
  REAL ( KIND= wp ) :: qsat_srf_old
  REAL ( KIND= wp ) :: evapo_pot
  INTEGER :: jsb_hori

!$acc data &
!$acc present(oro_stddev,t_srf_old,press_srf,drag_srf,q_air,wsr_max,w_skin &
!$acc ,w_snow_soil,w_snow_can,fract_snow_can,fract_water,fract_snow_soil)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( press_srf , 1 ) , 1
   IF ( wsr_max ( jsb_hori ) > epsilon ( 1._wp ) ) THEN
    fract_snow_can ( jsb_hori ) = min ( 1._wp , w_snow_can ( jsb_hori ) /&
     wsr_max ( jsb_hori ) )
    fract_water ( jsb_hori ) = min ( 1._wp , w_skin ( jsb_hori ) / wsr_max (&
     jsb_hori ) )
   ELSE
    fract_snow_can ( jsb_hori ) = 0._wp
    fract_water ( jsb_hori ) = 0._wp
   END IF
   fract_snow_soil ( jsb_hori ) = get_snow_fract_noforest ( w_snow_soil (&
    jsb_hori ) , oro_stddev ( jsb_hori ) )
   fract_snow_soil ( jsb_hori ) = merge ( fract_snow_can ( jsb_hori ) ,&
    fract_snow_soil ( jsb_hori ) , fract_snow_soil ( jsb_hori ) < epsilon (&
    1._wp ) .AND. fract_snow_can ( jsb_hori ) > epsilon ( 1._wp ) )
   qsat_srf_old = qsat ( t_srf_old ( jsb_hori ) , press_srf ( jsb_hori ) )
   evapo_pot = ( - 1._wp * heat_transfer_coef ( drag_srf ( jsb_hori ) ,&
    steplen ) * ( q_air ( jsb_hori ) - qsat_srf_old ) )
   IF ( fract_snow_soil ( jsb_hori ) > 0._wp ) THEN
    fract_snow_soil ( jsb_hori ) = fract_snow_soil ( jsb_hori ) / max ( 1._wp&
     , fract_snow_soil ( jsb_hori ) * evapo_pot * dtime / ( rhoh2o * (&
     w_snow_soil ( jsb_hori ) + w_snow_can ( jsb_hori ) ) ) )
   END IF
   IF ( fract_water ( jsb_hori ) > 0._wp ) THEN
    fract_water ( jsb_hori ) = fract_water ( jsb_hori ) / max ( 1._wp , (&
     1._wp - fract_snow_soil ( jsb_hori ) ) * evapo_pot * dtime / ( rhoh2o *&
     max ( epsilon ( 1._wp ) , w_skin ( jsb_hori ) ) ) )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_wskin_fractions_veg

 SUBROUTINE calc_wskin_fractions_bare ( dtime , steplen , wsr_max , oro_stddev&
  , t_srf , press_srf , drag_srf , q_air , w_skin , w_snow_soil , fract_water&
  , fract_snow_soil )
  USE mo_phy_schemes , ONLY: heat_transfer_coef , qsat
  USE mo_jsb_physical_constants , ONLY: rhoh2o

  REAL ( KIND= wp ) , INTENT(IN) :: dtime
  REAL ( KIND= wp ) , INTENT(IN) :: steplen
  REAL ( KIND= wp ) , INTENT(IN) :: oro_stddev ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: press_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: drag_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: q_air ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: wsr_max ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: w_skin ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: w_snow_soil ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: fract_water ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: fract_snow_soil ( : )
  REAL ( KIND= wp ) :: t_srf_old
  REAL ( KIND= wp ) :: qsat_srf_old
  REAL ( KIND= wp ) :: evapo_pot
  INTEGER :: jsb_hori

!$acc data &
!$acc present(oro_stddev,t_srf,press_srf,drag_srf,q_air,wsr_max,w_skin &
!$acc ,w_snow_soil,fract_water,fract_snow_soil)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( press_srf , 1 ) , 1
   IF ( wsr_max ( jsb_hori ) > epsilon ( 1._wp ) ) THEN
    fract_water ( jsb_hori ) = min ( 1._wp , w_skin ( jsb_hori ) / wsr_max (&
     jsb_hori ) )
   ELSE
    fract_water ( jsb_hori ) = 0._wp
   END IF
   fract_snow_soil ( jsb_hori ) = get_snow_fract_noforest ( w_snow_soil (&
    jsb_hori ) , oro_stddev ( jsb_hori ) )
   t_srf_old = t_srf ( jsb_hori )
   qsat_srf_old = qsat ( t_srf_old , press_srf ( jsb_hori ) )
   evapo_pot = heat_transfer_coef ( drag_srf ( jsb_hori ) , steplen ) * (&
    qsat_srf_old - q_air ( jsb_hori ) )
   IF ( fract_snow_soil ( jsb_hori ) > 0._wp .AND. w_snow_soil ( jsb_hori ) >&
    epsilon ( 1._wp ) ) THEN
    fract_snow_soil ( jsb_hori ) = fract_snow_soil ( jsb_hori ) / max ( 1._wp&
     , fract_snow_soil ( jsb_hori ) * evapo_pot * dtime / ( rhoh2o * max (&
     epsilon ( 1._wp ) , w_snow_soil ( jsb_hori ) ) ) )
   END IF
   IF ( fract_water ( jsb_hori ) > 0._wp .AND. w_skin ( jsb_hori ) > epsilon (&
    1._wp ) ) THEN
    fract_water ( jsb_hori ) = fract_water ( jsb_hori ) / max ( 1._wp , (&
     1._wp - fract_snow_soil ( jsb_hori ) ) * evapo_pot * dtime / ( rhoh2o *&
     max ( epsilon ( 1._wp ) , w_skin ( jsb_hori ) ) ) )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_wskin_fractions_bare

 FUNCTION get_snow_fract_noforest ( snow , orodev )
  USE mo_hydro_constants , ONLY: wsn2fract_const , wsn2fract_eps ,&
   wsn2fract_sigfac
  REAL ( KIND= wp ) :: get_snow_fract_noforest
  REAL ( KIND= wp ) , INTENT(IN) :: snow
  REAL ( KIND= wp ) , INTENT(IN) :: orodev

!$acc routine seq
  get_snow_fract_noforest = wsn2fract_const * tanh ( snow * 100._wp ) * sqrt (&
   snow * 1000._wp / ( snow * 1000._wp + wsn2fract_eps + wsn2fract_sigfac *&
   orodev ) )
 END FUNCTION get_snow_fract_noforest

 FUNCTION get_canopy_cond_unstressed_simple ( lai , par ) RESULT(conductance)
  USE mo_hydro_constants , ONLY: k => conductance_k , a => conductance_a , b&
   => conductance_b , c => conductance_c
  REAL ( KIND= wp ) :: conductance
  REAL ( KIND= wp ) , INTENT(IN) :: lai
  REAL ( KIND= wp ) , INTENT(IN) :: par
  REAL ( KIND= wp ) :: d
  REAL ( KIND= wp ) :: zpar
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_process" //&
  ":get_canopy_cond_unstressed_simple"

!$acc routine seq
  zpar = max ( 1.e-10_wp , par )
  IF ( lai > epsilon ( 1._wp ) ) THEN
   d = ( a + b * c ) / ( c * zpar )
   conductance = ( log ( ( d * exp ( k * lai ) + 1._wp ) / ( d + 1._wp ) ) * b&
    / ( d * zpar ) - log ( ( d + exp ( ( - k * lai ) ) ) / ( d + 1._wp ) ) ) /&
    ( k * c )
  ELSE
   conductance = epsilon ( 1._wp )
  END IF
 END FUNCTION get_canopy_cond_unstressed_simple

 FUNCTION get_canopy_cond_stressed_simple ( cond_unstressed , water_stress ,&
  air_is_saturated ) RESULT(conductance)
  REAL ( KIND= wp ) :: conductance
  REAL ( KIND= wp ) , INTENT(IN) :: cond_unstressed
  REAL ( KIND= wp ) , INTENT(IN) :: water_stress
  LOGICAL , INTENT(IN) :: air_is_saturated
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_hydro_process" //&
  ":get_canopy_cond_stressed_simple"

!$acc routine seq
  IF ( air_is_saturated ) THEN
   conductance = epsilon ( 1._wp )
  ELSE
   conductance = cond_unstressed * water_stress
  END IF
 END FUNCTION get_canopy_cond_stressed_simple

 FUNCTION get_water_stress_factor ( w_soil , w_soil_max , w_soil_crit_fract ,&
  w_soil_wilt_fract ) RESULT(water_stress_factor)
  REAL ( KIND= wp ) , INTENT(IN) :: w_soil
  REAL ( KIND= wp ) , INTENT(IN) :: w_soil_max
  REAL ( KIND= wp ) , INTENT(IN) :: w_soil_crit_fract
  REAL ( KIND= wp ) , INTENT(IN) :: w_soil_wilt_fract
  REAL ( KIND= wp ) :: water_stress_factor
  REAL ( KIND= wp ) :: w_crit
  REAL ( KIND= wp ) :: w_wilt

!$acc routine seq
  w_crit = w_soil_max * w_soil_crit_fract
  w_wilt = w_soil_max * w_soil_wilt_fract
  IF ( w_crit - w_wilt > 0._wp ) THEN
   water_stress_factor = max ( 0._wp , min ( 1._wp , ( w_soil - w_wilt ) / (&
    w_crit - w_wilt ) ) )
  ELSE
   water_stress_factor = 0._wp
  END IF
 END FUNCTION get_water_stress_factor

END MODULE mo_hydro_process

