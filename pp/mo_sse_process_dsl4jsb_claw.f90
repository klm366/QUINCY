MODULE mo_sse_process
 USE mo_jsb_physical_constants , ONLY: rhoh2o
 USE mo_kind , ONLY: wp
 PUBLIC :: relative_humidity_soil
 PUBLIC :: init_soil_temperature
 PUBLIC :: calc_vol_heat_capacity
 PUBLIC :: calc_thermal_conductivity
 PUBLIC :: calc_soil_temperature
 PUBLIC :: calc_snow_temperature
 PUBLIC :: calc_soiltemp_old
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_sse_process"
 PUBLIC :: get_liquid_max

CONTAINS
 FUNCTION calc_vol_heat_capacity ( vol_heat_cap_soil_dry_mineral , porosity ,&
  fract_organic , fract_water , fract_ice , layer_depth , dz )&
  RESULT(vol_heat_cap)
  USE mo_sse_constants , ONLY: vol_hcap_water , vol_hcap_ice , vol_hcap_org ,&
   vol_hcap_air , vol_hcap_bedrock
  REAL ( KIND= wp ) , INTENT(IN) :: porosity
  REAL ( KIND= wp ) , INTENT(IN) :: vol_heat_cap_soil_dry_mineral
  REAL ( KIND= wp ) , INTENT(IN) :: fract_organic
  REAL ( KIND= wp ) , INTENT(IN) :: fract_water
  REAL ( KIND= wp ) , INTENT(IN) :: fract_ice
  REAL ( KIND= wp ) , INTENT(IN) :: layer_depth
  REAL ( KIND= wp ) , INTENT(IN) :: dz
  REAL ( KIND= wp ) :: vol_heat_cap

!$acc routine seq
  IF ( layer_depth > 0._wp ) THEN
   vol_heat_cap = ( 1._wp - porosity ) * ( ( 1._wp - fract_organic ) *&
    vol_heat_cap_soil_dry_mineral + fract_organic * vol_hcap_org ) +&
    fract_water * vol_hcap_water + fract_ice * vol_hcap_ice + max ( porosity -&
    fract_water - fract_ice , 0._wp ) * vol_hcap_air
   vol_heat_cap = ( vol_heat_cap * layer_depth + vol_hcap_bedrock * ( dz -&
    layer_depth ) ) / dz
  ELSE
   vol_heat_cap = vol_hcap_bedrock
  END IF
 END FUNCTION calc_vol_heat_capacity

 FUNCTION calc_thermal_conductivity ( heat_cond_soil_mineral , water , ice ,&
  porosity , vol_field_cap , fract_organic , layer_depth , dz )&
  RESULT(heat_cond)
  USE mo_jsb_physical_constants , ONLY: dens_soil
  USE mo_sse_constants , ONLY: hcond_org , hcond_dry_org , hcond_water ,&
   hcond_ice , hcond_bedrock
  USE mo_hydro_constants , ONLY: vol_porosity_org
  REAL ( KIND= wp ) , INTENT(IN) :: heat_cond_soil_mineral
  REAL ( KIND= wp ) , INTENT(IN) :: water
  REAL ( KIND= wp ) , INTENT(IN) :: ice
  REAL ( KIND= wp ) , INTENT(IN) :: porosity
  REAL ( KIND= wp ) , INTENT(IN) :: vol_field_cap
  REAL ( KIND= wp ) , INTENT(IN) :: fract_organic
  REAL ( KIND= wp ) , INTENT(IN) :: layer_depth
  REAL ( KIND= wp ) , INTENT(IN) :: dz
  REAL ( KIND= wp ) :: heat_cond
  REAL ( KIND= wp ) :: dens_bulk
  REAL ( KIND= wp ) :: porosity_mineral
  REAL ( KIND= wp ) :: hcond_dry_min
  REAL ( KIND= wp ) :: hcond_dry
  REAL ( KIND= wp ) :: vpor_sat
  REAL ( KIND= wp ) :: saturation
  REAL ( KIND= wp ) :: kersten
  REAL ( KIND= wp ) :: hcond_soil
  REAL ( KIND= wp ) :: hcond_sat
  REAL ( KIND= wp ) :: fract_water
  REAL ( KIND= wp ) :: fract_ice

!$acc routine seq
  IF ( layer_depth > 0._wp ) THEN
   porosity_mineral = ( porosity - fract_organic * vol_porosity_org ) / (&
    1._wp - fract_organic )
   dens_bulk = dens_soil * ( 1._wp - porosity_mineral )
   hcond_dry_min = ( 0.135_wp * dens_bulk + 64.7_wp ) / ( dens_soil - 0.947_wp&
    * dens_bulk )
   hcond_dry = ( 1._wp - fract_organic ) * hcond_dry_min + fract_organic *&
    hcond_dry_org
   IF ( water + ice <= epsilon ( 1._wp ) ) THEN
    heat_cond = hcond_dry
   ELSE
    hcond_soil = ( 1._wp - fract_organic ) * heat_cond_soil_mineral +&
     fract_organic * hcond_org
    saturation = ( water + ice ) / layer_depth / porosity
    IF ( ice > 0._wp ) THEN
     kersten = saturation
    ELSE
     IF ( saturation > 0.1_wp ) THEN
      kersten = log10 ( saturation ) + 1._wp
     ELSE
      kersten = 0._wp
     END IF
    END IF
    vpor_sat = porosity - vol_field_cap
    fract_water = water / ( 1._wp - vpor_sat ) / layer_depth
    fract_ice = ice / ( 1._wp - vpor_sat ) / layer_depth
    hcond_sat = hcond_soil ** ( 1._wp - fract_water - fract_ice ) *&
     hcond_water ** ( fract_water ) * hcond_ice ** ( fract_ice )
    heat_cond = ( hcond_sat - hcond_dry ) * kersten + hcond_dry
    heat_cond = ( heat_cond * layer_depth + hcond_bedrock * ( dz - layer_depth&
     ) ) / dz
    heat_cond = max ( heat_cond , hcond_dry )
   END IF
  ELSE
   heat_cond = hcond_bedrock
  END IF
 END FUNCTION calc_thermal_conductivity

 FUNCTION relative_humidity_soil ( ws , fc )
  USE mo_math_constants , ONLY: pi
  REAL ( KIND= wp ) :: relative_humidity_soil
  REAL ( KIND= wp ) , INTENT(IN) :: ws
  REAL ( KIND= wp ) , INTENT(IN) :: fc
  REAL ( KIND= wp ) :: moisture

!$acc routine seq
  moisture = min ( ws , fc )
  IF ( moisture > 0._wp ) THEN
   relative_humidity_soil = 0.5_wp * ( 1._wp - cos ( moisture * pi / fc ) )
  ELSE
   relative_humidity_soil = 0._wp
  END IF
 END FUNCTION relative_humidity_soil

 SUBROUTINE init_soil_temperature ( cmid , tclim , tsoil , tsurf )
  USE mo_jsb_time , ONLY: t_datetime , get_time_start , get_year ,&
   get_year_length , get_month_length , get_day_length , get_year_day
  USE mo_math_constants , ONLY: pi

  REAL ( KIND= wp ) , INTENT(IN) :: cmid ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: tclim ( : , : )
  REAL ( KIND= wp ) , INTENT(OUT) :: tsoil ( : , : )
  REAL ( KIND= wp ) , INTENT(OUT) :: tsurf ( : )
  INTEGER :: nsoil
  INTEGER :: ndim
  REAL ( KIND= wp ) :: zkap
  REAL ( KIND= wp ) :: zsqrt
  REAL ( KIND= wp ) , ALLOCATABLE :: zrange ( : )
  REAL ( KIND= wp ) , ALLOCATABLE :: zdmax ( : )
  INTEGER , ALLOCATABLE :: jmax ( : )
  TYPE ( t_datetime ) , POINTER :: inidate
  INTEGER :: iniyear
  INTEGER :: im
  INTEGER :: i
  INTEGER :: jl
  REAL ( KIND= wp ) :: year_len
  REAL ( KIND= wp ) :: day_in_year
  REAL ( KIND= wp ) :: month_len ( 1 : 12 )
  REAL ( KIND= wp ) :: month_mid ( 1 : 12 )
  REAL ( KIND= wp ) , PARAMETER :: cmid_offset = 0.07_wp
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_sse_process" //&
  ":init_soil_temperature"

  ndim = size ( tsoil , 1 )
  nsoil = size ( tsoil , 2 )
  ALLOCATE ( zrange ( ndim ) , jmax ( ndim ) , zdmax ( ndim ) )
  zkap = 7.5e-7_wp
  inidate => get_time_start ( )
  iniyear = get_year ( inidate )
  year_len = real ( get_year_length ( iniyear ) , wp )
  day_in_year = int ( get_year_day ( inidate ) )
  DO im = 1 , 12 , 1
   month_len ( im ) = real ( get_month_length ( iniyear , im ) , wp )
   month_mid ( im ) = month_len ( im ) * 0.5_wp
  END DO
  zsqrt = sqrt ( zkap * year_len * get_day_length ( ) / pi )
  jmax ( : ) = maxloc ( tclim ( : , : ) , dim = 2 )
  zrange ( : ) = maxval ( tclim ( : , : ) , dim = 2 ) - minval ( tclim ( : , :&
   ) , dim = 2 )
  zdmax ( : ) = 0._wp
  DO jl = 1 , ndim , 1
   DO im = 1 , jmax ( jl ) , 1
    zdmax ( jl ) = zdmax ( jl ) + month_len ( im )
   END DO
   zdmax ( jl ) = zdmax ( jl ) - month_mid ( jmax ( jl ) )
  END DO
  DO i = 1 , nsoil , 1
   tsoil ( : , i ) = sum ( tclim ( : , : ) , dim = 2 ) / 12._wp + 0.5_wp *&
    zrange ( : ) * exp ( ( - ( cmid ( i ) - cmid_offset ) / zsqrt ) ) * cos (&
    2._wp * pi * ( day_in_year - zdmax ( : ) ) / year_len - ( cmid ( i ) -&
    cmid_offset ) / zsqrt )
  END DO
  tsurf ( : ) = tsoil ( : , 1 )
  DEALLOCATE ( zrange , jmax , zdmax )
 END SUBROUTINE init_soil_temperature

 SUBROUTINE calc_soil_temperature ( nidx , nsoil , dz , mids , delta_time ,&
  lstart , l_freeze , l_supercool , t_srf_unfilt , vol_heat_cap , heat_cond ,&
  t_soil , t_soil_acoef , t_soil_bcoef , hcap_grnd , grnd_hflx , thaw_depth ,&
  ws_max , matrix_pot , bclapp , ws , ice , frozen , melt )
  USE mo_jsb_physical_constants , ONLY: tmelt , alf , rhoh2o , rhoi

  INTEGER , INTENT(IN) :: nidx
  INTEGER , INTENT(IN) :: nsoil
  REAL ( KIND= wp ) , INTENT(IN) :: delta_time
  LOGICAL , INTENT(IN) :: lstart
  LOGICAL , INTENT(IN) :: l_freeze
  LOGICAL , INTENT(IN) :: l_supercool
  REAL ( KIND= wp ) , INTENT(IN) :: dz ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: mids ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf_unfilt ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: vol_heat_cap ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: heat_cond ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: t_soil ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: t_soil_acoef ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: t_soil_bcoef ( : , : )
  REAL ( KIND= wp ) , INTENT(OUT) :: hcap_grnd ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: grnd_hflx ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: thaw_depth ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: ws_max ( : , : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: matrix_pot ( : , : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(IN) :: bclapp ( : , : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: ws ( : , : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(INOUT) :: ice ( : , : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(OUT) :: frozen ( : , : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(OUT) :: melt ( : , : )
  INTEGER :: jl
  INTEGER :: jk
  REAL ( KIND= wp ) :: z1 ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: zd1 ( 1 : nsoil )
  REAL ( KIND= wp ) :: zdz1 ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: zdz2 ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: heat_cap ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: liquid_max ( 1 : nidx , 1 : nsoil )
  INTEGER :: claw_induction_0
  INTEGER :: claw_induction_1
  INTEGER :: claw_induction_2
  INTEGER :: claw_induction_3
  INTEGER :: claw_induction_4
  INTEGER :: claw_induction_5
  INTEGER :: claw_induction_6
  INTEGER :: claw_induction_7
  INTEGER :: claw_induction_10
  INTEGER :: claw_induction_11
  INTEGER :: claw_induction_12
  INTEGER :: claw_induction_13
  INTEGER :: claw_induction_16
  INTEGER :: claw_induction_17
  INTEGER :: claw_induction_18

INTEGER, SAVE :: ser_counter = 1
!$acc data create(z1,zd1,zdz1,zdz2,heat_cap,liquid_max) &
!$acc present(dz,mids,t_srf_unfilt,vol_heat_cap,heat_cond,t_soil, &
!$acc t_soil_acoef,t_soil_bcoef,hcap_grnd,grnd_hflx,thaw_depth)
!$acc parallel loop gang vector
  DO jk = 1 , nsoil - 1 , 1
   zd1 ( jk ) = 1._wp / ( mids ( jk + 1 ) - mids ( jk ) )
  END DO
!$acc parallel loop gang vector collapse(2)
  DO jl = 1 , nidx , 1
   DO jk = 1 , nsoil , 1
    heat_cap ( jl , jk ) = dz ( jk ) * vol_heat_cap ( jl , jk )
   END DO
  END DO
  IF ( present ( ws ) ) THEN
!$acc data present(frozen(:,:),melt(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_0 = 1 , size ( frozen , 1 )
    DO claw_induction_1 = 1 , size ( frozen , 2 )
     frozen ( claw_induction_0 , claw_induction_1 ) = 0._wp
     melt ( claw_induction_0 , claw_induction_1 ) = 0._wp
    END DO
   END DO
!$acc end parallel
!$acc end data
  END IF
!$acc data present(thaw_depth(:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_18 = 1 , size ( thaw_depth , 1 )
   thaw_depth ( claw_induction_18 ) = 0._wp
  END DO
!$acc end parallel
!$acc end data
  IF ( ( .NOT. lstart ) ) THEN
!$acc parallel default(none) &
!$acc present(t_srf_unfilt,t_soil,t_soil_acoef,t_soil_bcoef)
!$acc loop gang vector
   DO jl = 1 , nidx , 1
    t_soil ( jl , 1 ) = t_srf_unfilt ( jl )
   END DO
!$acc loop seq
   DO jk = 1 , nsoil - 1 , 1
!$acc loop gang vector
    DO jl = 1 , nidx , 1
     t_soil ( jl , jk + 1 ) = t_soil_acoef ( jl , jk ) + t_soil_bcoef ( jl ,&
      jk ) * t_soil ( jl , jk )
    END DO
   END DO
!$acc end parallel
   IF ( l_freeze .AND. present ( ws ) ) THEN
    IF ( l_supercool ) THEN
!$acc data &
!$acc present(t_soil(:,:),matrix_pot(:,:),bclapp(:,:),ws_max(:,:) &
!$acc ,liquid_max(:,:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_12 = 1 , size ( liquid_max , 1 )
      DO claw_induction_13 = 1 , size ( liquid_max , 2 )
       liquid_max ( claw_induction_12 , claw_induction_13 ) = get_liquid_max (&
        t_soil ( claw_induction_12 , claw_induction_13 ) , ws_max (&
        claw_induction_12 , claw_induction_13 ) , matrix_pot (&
        claw_induction_12 , claw_induction_13 ) , bclapp ( claw_induction_12 ,&
        claw_induction_13 ) )
      END DO
     END DO
!$acc end parallel
!$acc end data
    ELSE
!$acc data present(liquid_max(:,:))
!$acc parallel
!$acc loop gang vector
     DO claw_induction_10 = 1 , size ( liquid_max , 1 )
      DO claw_induction_11 = 1 , size ( liquid_max , 2 )
       liquid_max ( claw_induction_10 , claw_induction_11 ) = 0._wp
      END DO
     END DO
!$acc end parallel
!$acc end data
    END IF
!$acc data present(frozen,melt,ws,ice,t_soil,heat_cap,liquid_max)
!$acc parallel loop gang vector collapse(2)
    DO jl = 1 , nidx , 1
     DO jk = 1 , nsoil , 1
      IF ( t_soil ( jl , jk ) < tmelt .AND. ws ( jl , jk ) > liquid_max ( jl ,&
       jk ) ) THEN
       frozen ( jl , jk ) = min ( ws ( jl , jk ) - liquid_max ( jl , jk ) ,&
        heat_cap ( jl , jk ) * ( tmelt - t_soil ( jl , jk ) ) / ( alf * rhoh2o&
        ) )
      ELSE
       IF ( t_soil ( jl , jk ) > tmelt .AND. ice ( jl , jk ) > 0._wp ) THEN
        melt ( jl , jk ) = min ( ice ( jl , jk ) , heat_cap ( jl , jk ) * (&
         t_soil ( jl , jk ) - tmelt ) / ( alf * rhoi ) )
       END IF
      END IF
     END DO
    END DO
!$acc end parallel
!$acc end data
!$acc data &
!$acc present(frozen(:,:),melt(:,:),ice(:,:),t_soil(:,:),heat_cap(:,:),ws(: &
!$acc ,:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_2 = 1 , size ( ws , 1 )
     DO claw_induction_3 = 1 , size ( ws , 2 )
      ws ( claw_induction_2 , claw_induction_3 ) = ws ( claw_induction_2 ,&
       claw_induction_3 ) - frozen ( claw_induction_2 , claw_induction_3 ) +&
       melt ( claw_induction_2 , claw_induction_3 )
      ice ( claw_induction_2 , claw_induction_3 ) = ice ( claw_induction_2 ,&
       claw_induction_3 ) + frozen ( claw_induction_2 , claw_induction_3 ) -&
       melt ( claw_induction_2 , claw_induction_3 )
      t_soil ( claw_induction_2 , claw_induction_3 ) = t_soil (&
       claw_induction_2 , claw_induction_3 ) + frozen ( claw_induction_2 ,&
       claw_induction_3 ) * alf * rhoh2o / heat_cap ( claw_induction_2 ,&
       claw_induction_3 ) - melt ( claw_induction_2 , claw_induction_3 ) * alf&
       * rhoi / heat_cap ( claw_induction_2 , claw_induction_3 )
     END DO
    END DO
!$acc end parallel
!$acc end data
!$acc data present(frozen(:,:),melt(:,:))
!$acc parallel
!$acc loop gang vector
    DO claw_induction_4 = 1 , size ( frozen , 1 )
     DO claw_induction_5 = 1 , size ( frozen , 2 )
      frozen ( claw_induction_4 , claw_induction_5 ) = frozen (&
       claw_induction_4 , claw_induction_5 ) * rhoh2o / delta_time
      melt ( claw_induction_4 , claw_induction_5 ) = melt ( claw_induction_4 ,&
       claw_induction_5 ) * rhoi / delta_time
     END DO
    END DO
!$acc end parallel
!$acc end data
   END IF
!$acc parallel
!$acc loop gang vector
   DO jl = 1 , nidx , 1
    thaw_depth ( jl ) = mids ( nsoil )
    IF ( t_soil ( jl , 1 ) <= tmelt ) THEN
     thaw_depth ( jl ) = 0._wp
    END IF
!$acc loop seq
    DO jk = 2 , nsoil , 1
     IF ( thaw_depth ( jl ) == mids ( nsoil ) .AND. t_soil ( jl , jk - 1 ) >&
      tmelt .AND. t_soil ( jl , jk ) <= tmelt ) THEN
      thaw_depth ( jl ) = mids ( jk - 1 ) + ( mids ( jk ) - mids ( jk - 1 ) )&
       * ( t_soil ( jl , jk - 1 ) - tmelt ) / ( t_soil ( jl , jk - 1 ) -&
       t_soil ( jl , jk ) )
     END IF
    END DO
   END DO
!$acc end parallel
  END IF
!$acc parallel default(none) present(zdz1,zdz2,zd1,heat_cap,heat_cond)
!$acc loop gang vector collapse(2)
  DO jk = 1 , nsoil , 1
   DO jl = 1 , nidx , 1
    zdz2 ( jl , jk ) = heat_cap ( jl , jk ) / delta_time
   END DO
  END DO
!$acc loop gang vector collapse(2)
  DO jk = 1 , nsoil - 1 , 1
   DO jl = 1 , nidx , 1
    zdz1 ( jl , jk ) = zd1 ( jk ) * heat_cond ( jl , jk )
   END DO
  END DO
!$acc end parallel
!$acc data present(z1(:,:))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_16 = 1 , size ( z1 , 1 )
   DO claw_induction_17 = 1 , size ( z1 , 2 )
    z1 ( claw_induction_16 , claw_induction_17 ) = 0._wp
   END DO
  END DO
!$acc end parallel
!$acc end data
!$acc parallel loop default(none) present(zdz1,zdz2,z1) gang vector
  DO jl = 1 , nidx , 1
   z1 ( jl , nsoil - 1 ) = zdz2 ( jl , nsoil ) + zdz1 ( jl , nsoil - 1 )
  END DO
!$acc data &
!$acc present(t_soil_acoef(:,nsoil-1),t_soil_bcoef(:,nsoil-1),zdz2(:,nsoil) &
!$acc ,t_soil(:,nsoil),zdz1(:,nsoil-1),z1(:,nsoil-1))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_6 = 1 , size ( t_soil_acoef , 1 )
   t_soil_acoef ( claw_induction_6 , nsoil - 1 ) = zdz2 ( claw_induction_6 ,&
    nsoil ) * t_soil ( claw_induction_6 , nsoil ) / z1 ( claw_induction_6 ,&
    nsoil - 1 )
   t_soil_bcoef ( claw_induction_6 , nsoil - 1 ) = zdz1 ( claw_induction_6 ,&
    nsoil - 1 ) / z1 ( claw_induction_6 , nsoil - 1 )
  END DO
!$acc end parallel
!$acc end data
!$acc data present(z1,zdz1, zdz2,t_soil_acoef,t_soil_bcoef,t_soil)
!$acc parallel
!$acc loop gang vector
  DO jl = 1 , nidx , 1
!$acc loop seq
   DO jk = nsoil - 1 , 2 , (-1)
    z1 ( jl , jk - 1 ) = 1._wp / ( zdz2 ( jl , jk ) + zdz1 ( jl , jk - 1 ) +&
     zdz1 ( jl , jk ) * ( 1._wp - t_soil_bcoef ( jl , jk ) ) )
    t_soil_acoef ( jl , jk - 1 ) = ( t_soil ( jl , jk ) * zdz2 ( jl , jk ) +&
     zdz1 ( jl , jk ) * t_soil_acoef ( jl , jk ) ) * z1 ( jl , jk - 1 )
    t_soil_bcoef ( jl , jk - 1 ) = zdz1 ( jl , jk - 1 ) * z1 ( jl , jk - 1 )
   END DO
  END DO
!$acc end parallel
!$acc end data
!$acc data &
!$acc present(t_soil(:,1),t_soil_bcoef(:,1),grnd_hflx(:),zdz1(:,1) &
!$acc ,hcap_grnd(:),t_soil_acoef(:,1),zdz2(:,1))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_7 = 1 , size ( grnd_hflx , 1 )
   grnd_hflx ( claw_induction_7 ) = zdz1 ( claw_induction_7 , 1 ) * (&
    t_soil_acoef ( claw_induction_7 , 1 ) + ( t_soil_bcoef ( claw_induction_7&
    , 1 ) - 1._wp ) * t_soil ( claw_induction_7 , 1 ) )
   hcap_grnd ( claw_induction_7 ) = zdz2 ( claw_induction_7 , 1 ) * delta_time&
    + delta_time * ( 1._wp - t_soil_bcoef ( claw_induction_7 , 1 ) ) * zdz1 (&
    claw_induction_7 , 1 )
  END DO
!$acc end parallel
!$acc end data
!$acc end data
ser_counter = ser_counter + 1
 END SUBROUTINE calc_soil_temperature

 SUBROUTINE calc_snow_temperature ( nidx , nsnow , dz_soil , delta_time ,&
  lstart , dz , itop_old , itop , t_srf_unfilt , vol_heat_cap , heat_cond ,&
  vol_heat_cap_soil , heat_cond_soil , t_soil , t_soil_acoef , t_soil_bcoef ,&
  t_snow , t_snow_acoef , t_snow_bcoef , t_soil_top , hcap_grnd , grnd_hflx )

  INTEGER , INTENT(IN) :: nidx
  INTEGER , INTENT(IN) :: nsnow
  REAL ( KIND= wp ) , INTENT(IN) :: delta_time
  LOGICAL , INTENT(IN) :: lstart
  REAL ( KIND= wp ) , INTENT(IN) :: dz ( : , : )
  REAL ( KIND= wp ) , INTENT(IN) :: dz_soil ( 1 : 2 )
  INTEGER , INTENT(IN) :: itop ( : )
  INTEGER , INTENT(IN) :: itop_old ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf_unfilt ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: vol_heat_cap ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: heat_cond ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: vol_heat_cap_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: heat_cond_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_soil_acoef ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_soil_bcoef ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: t_snow ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: t_snow_acoef ( : , : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: t_snow_bcoef ( : , : )
  REAL ( KIND= wp ) , INTENT(OUT) :: t_soil_top ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: hcap_grnd ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: grnd_hflx ( : )
  INTEGER :: jl
  INTEGER :: jk
  INTEGER :: itop_new
  REAL ( KIND= wp ) :: zmid ( 1 : nidx , 1 : nsnow + 2 )
  REAL ( KIND= wp ) :: z1 ( 1 : nidx )
  REAL ( KIND= wp ) :: zd1 ( 1 : nidx , 1 : nsnow + 1 )
  REAL ( KIND= wp ) :: zdz1 ( 1 : nidx , 1 : nsnow + 1 )
  REAL ( KIND= wp ) :: zdz2 ( 1 : nidx , 1 : nsnow + 1 )
  INTEGER :: claw_induction_8
  INTEGER :: claw_induction_9
  INTEGER :: claw_induction_14
  INTEGER :: claw_induction_15

!$acc data create(zmid,z1,zd1,zdz1,zdz2)
!$acc data &
!$acc present(dz(:,1),grnd_hflx(:),t_soil_top(:),hcap_grnd(:),t_srf_unfilt(:) &
!$acc ,zmid(:,1))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_8 = 1 , size ( grnd_hflx , 1 )
   grnd_hflx ( claw_induction_8 ) = 0._wp
   hcap_grnd ( claw_induction_8 ) = 0._wp
   t_soil_top ( claw_induction_8 ) = t_srf_unfilt ( claw_induction_8 )
   zmid ( claw_induction_8 , 1 ) = dz ( claw_induction_8 , 1 ) / 2._wp
  END DO
!$acc end parallel
!$acc end data
!$acc data present(zmid,dz,zd1)
!$acc parallel
!$acc loop gang vector
  DO jl = 1 , nidx , 1
!$acc loop seq
   DO jk = 2 , nsnow , 1
    zmid ( jl , jk ) = zmid ( jl , jk - 1 ) + 0.5_wp * ( dz ( jl , jk - 1 ) +&
     dz ( jl , jk ) )
    IF ( zmid ( jl , jk ) > zmid ( jl , jk - 1 ) ) THEN
     zd1 ( jl , jk - 1 ) = 1._wp / ( zmid ( jl , jk ) - zmid ( jl , jk - 1 ) )
    END IF
   END DO
  END DO
!$acc end parallel
!$acc end data
!$acc data &
!$acc present(zmid(:,nsnow),dz_soil(1),zmid(:,nsnow+1),dz(:,nsnow),zmid(: &
!$acc ,nsnow+2),zd1(:,nsnow),dz_soil(2),zd1(:,nsnow+1))
!$acc parallel
!$acc loop gang vector
  DO claw_induction_9 = 1 , size ( zmid , 1 )
   zmid ( claw_induction_9 , nsnow + 1 ) = zmid ( claw_induction_9 , nsnow ) +&
    0.5_wp * ( dz ( claw_induction_9 , nsnow ) + dz_soil ( 1 ) )
   zmid ( claw_induction_9 , nsnow + 2 ) = zmid ( claw_induction_9 , nsnow + 1&
    ) + 0.5_wp * ( dz_soil ( 1 ) + dz_soil ( 2 ) )
   zd1 ( claw_induction_9 , nsnow ) = 1._wp / ( zmid ( claw_induction_9 ,&
    nsnow + 1 ) - zmid ( claw_induction_9 , nsnow ) )
   zd1 ( claw_induction_9 , nsnow + 1 ) = 1._wp / ( zmid ( claw_induction_9 ,&
    nsnow + 2 ) - zmid ( claw_induction_9 , nsnow + 1 ) )
  END DO
!$acc end parallel
!$acc end data
!$acc data present(itop,t_snow)
!$acc parallel
!$acc loop gang vector private(itop_new)
  DO jl = 1 , nidx , 1
   itop_new = itop ( jl )
!$acc loop seq
   DO jk = 1 , itop_new - 1 , 1
    t_snow ( jl , jk ) = 0._wp
   END DO
  END DO
!$acc end parallel
!$acc end data
  IF ( ( .NOT. lstart ) ) THEN
!$acc data present(t_snow(:,:))
!$acc parallel
!$acc loop gang vector
   DO claw_induction_14 = 1 , size ( t_snow , 1 )
    DO claw_induction_15 = 1 , size ( t_snow , 2 )
     t_snow ( claw_induction_14 , claw_induction_15 ) = 0._wp
    END DO
   END DO
!$acc end parallel
!$acc end data
!$acc data &
!$acc present(t_snow,itop,itop_old,t_srf_unfilt,t_snow_acoef,t_snow_bcoef &
!$acc ,t_soil_top)
!$acc parallel loop gang vector collapse(2)
   DO jl = 1 , nidx , 1
    DO jk = 1 , nsnow , 1
     IF ( itop ( jl ) <= jk .AND. jk <= itop_old ( jl ) ) THEN
      t_snow ( jl , jk ) = t_srf_unfilt ( jl )
     ELSE
      IF ( itop ( jl ) > itop_old ( jl ) .AND. jk == itop ( jl ) ) THEN
       t_snow ( jl , jk ) = t_srf_unfilt ( jl )
      END IF
     END IF
    END DO
   END DO
!$acc parallel loop gang vector
   DO jl = 1 , nidx , 1
!$acc loop seq
    DO jk = 2 , nsnow , 1
     IF ( jk > max ( itop ( jl ) , itop_old ( jl ) ) ) THEN
      t_snow ( jl , jk ) = t_snow_acoef ( jl , jk - 1 ) + t_snow_bcoef ( jl ,&
       jk - 1 ) * t_snow ( jl , jk - 1 )
     END IF
    END DO
   END DO
!$acc parallel loop gang vector
   DO jl = 1 , nidx , 1
    IF ( max ( itop ( jl ) , itop_old ( jl ) ) <= nsnow ) THEN
     t_soil_top ( jl ) = t_snow_acoef ( jl , nsnow ) + t_snow_bcoef ( jl ,&
      nsnow ) * t_snow ( jl , nsnow )
    ELSE
     t_soil_top ( jl ) = t_srf_unfilt ( jl )
    END IF
   END DO
!$acc end data
  END IF
!$acc data present(itop,t_snow_acoef,t_snow_bcoef)
!$acc parallel
!$acc loop gang vector private(itop_new)
  DO jl = 1 , nidx , 1
   itop_new = itop ( jl )
!$acc loop seq
   DO jk = 1 , itop_new - 1 , 1
    t_snow_acoef ( jl , jk ) = 0._wp
    t_snow_bcoef ( jl , jk ) = 0._wp
   END DO
  END DO
!$acc end parallel
!$acc end data
!$acc parallel loop gang vector collapse(2)
  DO jl = 1 , nidx , 1
   DO jk = 1 , nsnow , 1
    IF ( jk >= itop ( jl ) ) THEN
     zdz2 ( jl , jk ) = vol_heat_cap ( jl ) * dz ( jl , jk ) / delta_time
     zdz1 ( jl , jk ) = zd1 ( jl , jk ) * heat_cond ( jl )
    END IF
   END DO
  END DO
!$acc data &
!$acc present(zdz1,zdz2,vol_heat_cap_soil,dz_soil,zd1,heat_cond_soil,z1 &
!$acc ,t_soil_bcoef,t_snow_acoef,t_snow_bcoef,t_soil,t_soil_acoef)
!$acc parallel loop gang vector
  DO jl = 1 , nidx , 1
   zdz2 ( jl , nsnow + 1 ) = vol_heat_cap_soil ( jl ) * dz_soil ( 1 ) /&
    delta_time
   zdz1 ( jl , nsnow + 1 ) = zd1 ( jl , nsnow + 1 ) * heat_cond_soil ( jl )
   IF ( itop ( jl ) <= nsnow ) THEN
    z1 ( jl ) = 1._wp / ( zdz2 ( jl , nsnow + 1 ) + zdz1 ( jl , nsnow ) + zdz1&
     ( jl , nsnow + 1 ) * ( 1._wp - t_soil_bcoef ( jl ) ) )
    t_snow_acoef ( jl , nsnow ) = ( t_soil ( jl ) * zdz2 ( jl , nsnow + 1 ) +&
     zdz1 ( jl , nsnow + 1 ) * t_soil_acoef ( jl ) ) * z1 ( jl )
    t_snow_bcoef ( jl , nsnow ) = zdz1 ( jl , nsnow ) * z1 ( jl )
   ELSE
    t_snow_acoef ( jl , nsnow ) = 0._wp
    t_snow_bcoef ( jl , nsnow ) = 0._wp
   END IF
  END DO
!$acc end data
!$acc data present(itop,z1,zdz1,zdz2,t_snow_acoef,t_snow_bcoef,t_snow)
!$acc parallel loop gang vector
  DO jl = 1 , nidx , 1
!$acc loop seq
   DO jk = nsnow , 2 , (-1)
    IF ( jk > itop ( jl ) ) THEN
     z1 ( jl ) = 1._wp / ( zdz2 ( jl , jk ) + zdz1 ( jl , jk - 1 ) + zdz1 ( jl&
      , jk ) * ( 1._wp - t_snow_bcoef ( jl , jk ) ) )
     t_snow_acoef ( jl , jk - 1 ) = ( t_snow ( jl , jk ) * zdz2 ( jl , jk ) +&
      zdz1 ( jl , jk ) * t_snow_acoef ( jl , jk ) ) * z1 ( jl )
     t_snow_bcoef ( jl , jk - 1 ) = zdz1 ( jl , jk - 1 ) * z1 ( jl )
    END IF
   END DO
  END DO
!$acc end data
!$acc data &
!$acc present(itop,grnd_hflx,hcap_grnd,zdz1,zdz2,t_snow_acoef,t_snow_bcoef &
!$acc ,t_snow)
!$acc parallel loop gang vector collapse(2)
  DO jl = 1 , nidx , 1
   DO jk = 1 , nsnow , 1
    IF ( itop ( jl ) == jk ) THEN
     grnd_hflx ( jl ) = zdz1 ( jl , jk ) * ( t_snow_acoef ( jl , jk ) + (&
      t_snow_bcoef ( jl , jk ) - 1._wp ) * t_snow ( jl , jk ) )
     hcap_grnd ( jl ) = zdz2 ( jl , jk ) * delta_time + delta_time * ( 1._wp -&
      t_snow_bcoef ( jl , jk ) ) * zdz1 ( jl , jk )
    END IF
   END DO
  END DO
!$acc end data
!$acc end data
 END SUBROUTINE calc_snow_temperature

 SUBROUTINE calc_soiltemp_old ( nidx , nsoil , cdel , cmid , delta_time ,&
  lstart , pts , psn , psocond , prgcgn , pgrndc , pgrndd , ptsoil , pgrndcapc&
  , pgrndhflx , ldglac )

  INTEGER , INTENT(IN) :: nidx
  INTEGER , INTENT(IN) :: nsoil
  REAL ( KIND= wp ) , INTENT(IN) :: cdel ( 1 : nsoil )
  REAL ( KIND= wp ) , INTENT(IN) :: cmid ( 1 : nsoil )
  REAL ( KIND= wp ) , INTENT(IN) :: delta_time
  LOGICAL , INTENT(IN) :: lstart
  REAL ( KIND= wp ) , INTENT(IN) :: pts ( 1 : nidx )
  REAL ( KIND= wp ) , INTENT(IN) :: psn ( 1 : nidx )
  REAL ( KIND= wp ) , INTENT(IN) :: psocond ( 1 : nidx )
  REAL ( KIND= wp ) , INTENT(IN) :: prgcgn ( 1 : nidx )
  REAL ( KIND= wp ) , INTENT(INOUT) :: pgrndc ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) , INTENT(INOUT) :: pgrndd ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) , INTENT(INOUT) :: ptsoil ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) , INTENT(OUT) :: pgrndcapc ( 1 : nidx )
  REAL ( KIND= wp ) , INTENT(OUT) :: pgrndhflx ( 1 : nidx )
  LOGICAL , INTENT(IN) :: ldglac ( 1 : nidx )
  INTEGER :: jk
  INTEGER :: i
  REAL ( KIND= wp ) :: zso_cond ( 1 : nidx )
  REAL ( KIND= wp ) :: zso_capa ( 1 : nidx )
  REAL ( KIND= wp ) :: z1 ( 1 : nidx )
  REAL ( KIND= wp ) :: zd1 ( 1 : nsoil )
  REAL ( KIND= wp ) :: zdz1 ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: zdz2 ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: zkappa ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: zcapa ( 1 : nidx , 1 : nsoil )
  REAL ( KIND= wp ) :: zsnow_h ( 1 : nidx )
  REAL ( KIND= wp ) :: zx1 ( 1 : nidx )
  REAL ( KIND= wp ) :: zx2 ( 1 : nidx )
  REAL ( KIND= wp ) :: zrici
  REAL ( KIND= wp ) :: zdifiz
  REAL ( KIND= wp ) :: zsn_cond
  REAL ( KIND= wp ) :: zsn_dens
  REAL ( KIND= wp ) :: zsn_capa

  zrici = 2.09e+06_wp
  zdifiz = 12.e-07_wp
  zsn_cond = 0.31_wp
  zsn_dens = 330.0_wp
  zsn_capa = 634500.0_wp
!$acc data present(cdel,cmid,pts,psn,psocond,prgcgn,pgrndc,pgrndd,ptsoil,pgrndcapc,pgrndhflx,ldglac) &
!$acc create(zso_cond,zso_capa,z1,zd1,zdz1,zdz2,zkappa,zcapa,zsnow_h,zx1,zx2)
!$acc parallel loop gang vector
  DO jk = 1 , nsoil - 1 , 1
   zd1 ( jk ) = 1._wp / ( cmid ( jk + 1 ) - cmid ( jk ) )
  END DO
!$acc parallel loop gang vector
  DO i = 1 , nidx , 1
   IF ( ldglac ( i ) ) THEN
    zso_capa ( i ) = zrici
    zso_cond ( i ) = zso_capa ( i ) * zdifiz
   ELSE
    zso_capa ( i ) = prgcgn ( i )
    zso_cond ( i ) = psocond ( i )
   END IF
  END DO
!$acc parallel loop gang vector collapse(2)
  DO i = 1 , nidx , 1
   DO jk = 1 , nsoil , 1
    zkappa ( i , jk ) = zso_cond ( i )
    zcapa ( i , jk ) = zso_capa ( i )
   END DO
  END DO
  IF ( ( .NOT. lstart ) ) THEN
!$acc parallel loop gang vector
   DO i = 1 , nidx , 1
    ptsoil ( i , 1 ) = pts ( i )
   END DO
!$acc parallel loop gang vector collapse(2)
   DO i = 1 , nidx , 1
    DO jk = 1 , nsoil - 1 , 1
     ptsoil ( i , jk + 1 ) = pgrndc ( i , jk ) + pgrndd ( i , jk ) * ptsoil (&
      i , jk )
    END DO
   END DO
  END IF
!$acc parallel loop gang vector
  DO i = 1 , nidx , 1
   zsnow_h ( i ) = psn ( i ) * rhoh2o / zsn_dens
   IF ( zsnow_h ( i ) > cmid ( 2 ) ) THEN
    zcapa ( i , 1 ) = zsn_capa
    zkappa ( i , 1 ) = zsn_cond
   ELSE
    IF ( zsnow_h ( i ) > 0.0_wp .AND. zsnow_h ( i ) <= cmid ( 2 ) ) THEN
     zx1 ( i ) = zsnow_h ( i ) / cmid ( 2 )
     zx2 ( i ) = ( cmid ( 2 ) - zsnow_h ( i ) ) / cmid ( 2 )
     zcapa ( i , 1 ) = zx1 ( i ) * zsn_capa + zx2 ( i ) * zso_capa ( i )
     zkappa ( i , 1 ) = 1._wp / ( zx1 ( i ) / zsn_cond + zx2 ( i ) / zso_cond&
      ( i ) )
    ELSE
     zcapa ( i , 1 ) = zso_capa ( i )
     zkappa ( i , 1 ) = zso_cond ( i )
    END IF
   END IF
  END DO
!$acc parallel
!$acc loop gang vector
  DO i = 1 , nidx , 1
!$acc loop seq
   DO jk = 2 , nsoil - 2 , 1
    IF ( zsnow_h ( i ) > cmid ( jk + 1 ) ) THEN
     zcapa ( i , jk ) = zsn_capa
     zkappa ( i , jk ) = zsn_cond
    ELSE
     IF ( zsnow_h ( i ) > cmid ( jk ) .AND. zsnow_h ( i ) <= cmid ( jk + 1 ) )&
      THEN
      zx1 ( i ) = ( zsnow_h ( i ) - cmid ( jk ) ) * zd1 ( jk )
      zx2 ( i ) = ( cmid ( jk + 1 ) - zsnow_h ( i ) ) * zd1 ( jk )
      zcapa ( i , jk ) = zx1 ( i ) * zsn_capa + zx2 ( i ) * zso_capa ( i )
      zkappa ( i , jk ) = 1._wp / ( zx1 ( i ) / zsn_cond + zx2 ( i ) /&
       zso_cond ( i ) )
     ELSE
      zcapa ( i , jk ) = zso_capa ( i )
      zkappa ( i , jk ) = zso_cond ( i )
     END IF
    END IF
   END DO
  END DO
!$acc end parallel
!$acc parallel loop gang vector collapse(2)
  DO i = 1 , nidx , 1
   DO jk = 1 , nsoil , 1
    zdz2 ( i , jk ) = zcapa ( i , jk ) * cdel ( jk ) / delta_time
   END DO
  END DO
!$acc parallel loop gang vector collapse(2)
  DO i = 1 , nidx , 1
   DO jk = 1 , nsoil - 1 , 1
    zdz1 ( i , jk ) = zd1 ( jk ) * zkappa ( i , jk )
   END DO
  END DO
!$acc parallel loop gang vector
  DO i = 1 , nidx , 1
   z1 ( i ) = zdz2 ( i , nsoil ) + zdz1 ( i , nsoil - 1 )
   pgrndc ( i , nsoil - 1 ) = zdz2 ( i , nsoil ) * ptsoil ( i , nsoil ) / z1 (&
    i )
   pgrndd ( i , nsoil - 1 ) = zdz1 ( i , nsoil - 1 ) / z1 ( i )
  END DO
!$acc parallel
!$acc loop gang vector
  DO i = 1 , nidx , 1
!$acc loop seq
   DO jk = nsoil - 1 , 2 , (-1)
    z1 ( i ) = 1._wp / ( zdz2 ( i , jk ) + zdz1 ( i , jk - 1 ) + zdz1 ( i , jk&
     ) * ( 1._wp - pgrndd ( i , jk ) ) )
    pgrndc ( i , jk - 1 ) = ( ptsoil ( i , jk ) * zdz2 ( i , jk ) + zdz1 ( i ,&
     jk ) * pgrndc ( i , jk ) ) * z1 ( i )
    pgrndd ( i , jk - 1 ) = zdz1 ( i , jk - 1 ) * z1 ( i )
   END DO
  END DO
!$acc end parallel
!$acc parallel loop gang vector
  DO i = 1 , nidx , 1
   pgrndhflx ( i ) = zdz1 ( i , 1 ) * ( pgrndc ( i , 1 ) + ( pgrndd ( i , 1 )&
    - 1._wp ) * ptsoil ( i , 1 ) )
   pgrndcapc ( i ) = zdz2 ( i , 1 ) * delta_time + delta_time * ( 1._wp -&
    pgrndd ( i , 1 ) ) * zdz1 ( i , 1 )
  END DO
!$acc end data
 END SUBROUTINE calc_soiltemp_old

 FUNCTION get_liquid_max ( temp , w_max , matrix_pot , bclapp )&
  RESULT(liquid_max)
  USE mo_jsb_physical_constants , ONLY: tmelt , alf , grav
  REAL ( KIND= wp ) , INTENT(IN) :: temp
  REAL ( KIND= wp ) , INTENT(IN) :: w_max
  REAL ( KIND= wp ) , INTENT(IN) :: matrix_pot
  REAL ( KIND= wp ) , INTENT(IN) :: bclapp
  REAL ( KIND= wp ) :: liquid_max

!$acc routine seq
  IF ( temp < tmelt ) THEN
   liquid_max = w_max * ( alf * ( tmelt - temp ) / max ( 1.e-10_wp , grav *&
    temp * matrix_pot ) ) ** ( ( - 1._wp / max ( 1._wp , bclapp ) ) )
  ELSE
   liquid_max = 0._wp
  END IF
 END FUNCTION get_liquid_max

END MODULE mo_sse_process

