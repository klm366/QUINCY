MODULE mo_pheno_process
 USE mo_kind , ONLY: wp
 PRIVATE :: get_letitdie
 PRIVATE :: get_letitgrow
 PRIVATE :: get_shedrate_rg
 PUBLIC :: calc_summergreen_phenology
 PUBLIC :: calc_evergreen_phenology
 PUBLIC :: calc_crop_phenology
 PUBLIC :: calc_raingreen_phenology
 PUBLIC :: calc_grass_phenology
 PUBLIC :: get_foliage_projected_cover
 PUBLIC :: derive_maxlai_according_to_allometric_relationships
 PUBLIC :: track_max_green_pool
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_pheno_process"

CONTAINS
 SUBROUTINE calc_summergreen_phenology ( day_of_year , is_newday ,&
  timestep_in_days , lat , previous_day_temp_mean , pseudo_soil_temp ,&
  laimax_corrected , days_since_growth_begin_sg , chill_days_sg , heat_sum_sg&
  , growth_phase_sg , lai )
  USE mo_pheno_constants , ONLY: phenoparam

  INTEGER , INTENT(IN) :: day_of_year
  LOGICAL , INTENT(IN) :: is_newday
  REAL ( KIND= wp ) , INTENT(IN) :: timestep_in_days
  REAL ( KIND= wp ) , INTENT(IN) :: lat ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: previous_day_temp_mean ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: laimax_corrected ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: days_since_growth_begin_sg ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: chill_days_sg ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: heat_sum_sg ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: growth_phase_sg ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: lai ( : )
  REAL ( KIND= wp ) :: xtmp
  REAL ( KIND= wp ) :: grrate
  INTEGER :: jsb_hori

!$acc data &
!$acc present(lat,previous_day_temp_mean,laimax_corrected,pseudo_soil_temp &
!$acc ,days_since_growth_begin_sg,chill_days_sg,heat_sum_sg,growth_phase_sg &
!$acc ,lai)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( growth_phase_sg , 1 ) , 1
   IF ( growth_phase_sg ( jsb_hori ) < ( - 0.5_wp ) .AND. heat_sum_sg (&
    jsb_hori ) < ( - 1._wp ) .AND. ( lat ( jsb_hori ) >= 0._wp .AND.&
    day_of_year == 1 .OR. lat ( jsb_hori ) < 0._wp .AND. day_of_year == 183 )&
    ) THEN
    heat_sum_sg ( jsb_hori ) = 0._wp
   END IF
   IF ( is_newday ) THEN
    IF ( days_since_growth_begin_sg ( jsb_hori ) > 0.0_wp .AND.&
     days_since_growth_begin_sg ( jsb_hori ) < 365.0_wp ) THEN
     days_since_growth_begin_sg ( jsb_hori ) = days_since_growth_begin_sg (&
      jsb_hori ) + 1.0_wp
    END IF
    IF ( heat_sum_sg ( jsb_hori ) > ( - 1._wp ) ) THEN
     IF ( pseudo_soil_temp ( jsb_hori ) < phenoparam % sg % alternation_temp )&
      THEN
      chill_days_sg ( jsb_hori ) = chill_days_sg ( jsb_hori ) + 1._wp
      chill_days_sg ( jsb_hori ) = min ( real ( chill_days_sg ( jsb_hori ) ) ,&
       real ( phenoparam % eg_sg % max_chill_days ) )
     ELSE
      heat_sum_sg ( jsb_hori ) = heat_sum_sg ( jsb_hori ) + pseudo_soil_temp (&
       jsb_hori ) - phenoparam % sg % alternation_temp
     END IF
    END IF
    IF ( growth_phase_sg ( jsb_hori ) < ( - 0.5_wp ) ) THEN
     IF ( heat_sum_sg ( jsb_hori ) > ( - 1._wp ) .AND. heat_sum_sg ( jsb_hori&
      ) - phenoparam % sg % heat_sum_min > phenoparam % sg % heat_sum_range *&
      exp ( ( - chill_days_sg ( jsb_hori ) / phenoparam % sg %&
      chill_decay_const ) ) ) THEN
      growth_phase_sg ( jsb_hori ) = 0._wp
      days_since_growth_begin_sg ( jsb_hori ) = 1.0_wp
      heat_sum_sg ( jsb_hori ) = ( - 99._wp )
      chill_days_sg ( jsb_hori ) = 0._wp
     END IF
    ELSE
     IF ( growth_phase_sg ( jsb_hori ) > 0.5_wp ) THEN
      IF ( previous_day_temp_mean ( jsb_hori ) < pseudo_soil_temp ( jsb_hori )&
       .AND. pseudo_soil_temp ( jsb_hori ) < phenoparam % sg %&
       autumn_event_temp .AND. days_since_growth_begin_sg ( jsb_hori ) >&
       phenoparam % sg % minlength_gvphase ) THEN
       growth_phase_sg ( jsb_hori ) = ( - 1._wp )
      ELSE
       IF ( days_since_growth_begin_sg ( jsb_hori ) > phenoparam % sg %&
        maxlength_gvphase ) THEN
        growth_phase_sg ( jsb_hori ) = ( - 1._wp )
       END IF
      END IF
     ELSE
      IF ( growth_phase_sg ( jsb_hori ) > ( - 0.5_wp ) .AND. growth_phase_sg (&
       jsb_hori ) < 0.5_wp ) THEN
       IF ( days_since_growth_begin_sg ( jsb_hori ) > phenoparam % sg %&
        growthphaselength ) THEN
        growth_phase_sg ( jsb_hori ) = 1.0_wp
       END IF
      END IF
     END IF
    END IF
   END IF
   IF ( growth_phase_sg ( jsb_hori ) < ( - 0.5_wp ) ) THEN
    xtmp = phenoparam % sg % shedrate_rest
    lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
     xtmp )
   ELSE
    IF ( growth_phase_sg ( jsb_hori ) > 0.5_wp ) THEN
     xtmp = phenoparam % sg % shedrate_veget
     lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
      xtmp )
    ELSE
     IF ( growth_phase_sg ( jsb_hori ) > ( - 0.5_wp ) .AND. growth_phase_sg (&
      jsb_hori ) < 0.5_wp ) THEN
      IF ( lai ( jsb_hori ) < phenoparam % all % laiseed ) THEN
       lai ( jsb_hori ) = min ( laimax_corrected ( jsb_hori ) - epsilon (&
        1.0_wp ) , phenoparam % all % laiseed )
      END IF
      IF ( lai ( jsb_hori ) > phenoparam % all % lai_negligible ) THEN
       grrate = phenoparam % sg % growthrate
       lai ( jsb_hori ) = get_letitgrow ( timestep_in_days , lai ( jsb_hori )&
        , grrate , 0.0_wp , laimax_corrected ( jsb_hori ) )
      END IF
     END IF
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_summergreen_phenology

 SUBROUTINE calc_evergreen_phenology ( day_of_year , is_newday ,&
  timestep_in_days , lat , pseudo_soil_temp , laimax_corrected ,&
  days_since_growth_begin_eg , chill_days_eg , heat_sum_eg , growth_phase_eg ,&
  lai )
  USE mo_pheno_constants , ONLY: phenoparam

  INTEGER , INTENT(IN) :: day_of_year
  LOGICAL , INTENT(IN) :: is_newday
  REAL ( KIND= wp ) , INTENT(IN) :: timestep_in_days
  REAL ( KIND= wp ) , INTENT(IN) :: lat ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: laimax_corrected ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: days_since_growth_begin_eg ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: chill_days_eg ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: heat_sum_eg ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: growth_phase_eg ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: lai ( : )
  REAL ( KIND= wp ) :: xtmp
  REAL ( KIND= wp ) :: grrate
  INTEGER :: jsb_hori

!$acc data &
!$acc present(lat,laimax_corrected,pseudo_soil_temp &
!$acc ,days_since_growth_begin_eg,chill_days_eg,heat_sum_eg,growth_phase_eg &
!$acc ,lai)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( laimax_corrected , 1 ) , 1
   IF ( growth_phase_eg ( jsb_hori ) > 0.5_wp .AND. heat_sum_eg ( jsb_hori ) <&
    ( - 1._wp ) .AND. ( lat ( jsb_hori ) >= 0._wp .AND. day_of_year == 1 .OR.&
    lat ( jsb_hori ) < 0._wp .AND. day_of_year == 183 ) ) THEN
    heat_sum_eg ( jsb_hori ) = 0._wp
   END IF
   IF ( is_newday ) THEN
    IF ( days_since_growth_begin_eg ( jsb_hori ) > 0.0_wp .AND.&
     days_since_growth_begin_eg ( jsb_hori ) < 365.0_wp ) THEN
     days_since_growth_begin_eg ( jsb_hori ) = days_since_growth_begin_eg (&
      jsb_hori ) + 1.0_wp
    END IF
    IF ( heat_sum_eg ( jsb_hori ) > ( - 1.0_wp ) ) THEN
     IF ( pseudo_soil_temp ( jsb_hori ) < phenoparam % eg % alternation_temp )&
      THEN
      chill_days_eg ( jsb_hori ) = chill_days_eg ( jsb_hori ) + 1._wp
      chill_days_eg ( jsb_hori ) = min ( real ( chill_days_eg ( jsb_hori ) ) ,&
       real ( phenoparam % eg_sg % max_chill_days ) )
     ELSE
      heat_sum_eg ( jsb_hori ) = heat_sum_eg ( jsb_hori ) + pseudo_soil_temp (&
       jsb_hori ) - phenoparam % eg % alternation_temp
     END IF
    END IF
    IF ( growth_phase_eg ( jsb_hori ) < 0.5_wp ) THEN
     IF ( days_since_growth_begin_eg ( jsb_hori ) > phenoparam % eg %&
      growthphaselength ) THEN
      growth_phase_eg ( jsb_hori ) = 1.0_wp
     END IF
    ELSE
     IF ( heat_sum_eg ( jsb_hori ) > ( - 1._wp ) .AND. heat_sum_eg ( jsb_hori&
      ) - phenoparam % eg % heat_sum_min > phenoparam % eg % heat_sum_range *&
      exp ( ( - chill_days_eg ( jsb_hori ) / phenoparam % eg %&
      chill_decay_const ) ) ) THEN
      growth_phase_eg ( jsb_hori ) = 0._wp
      days_since_growth_begin_eg ( jsb_hori ) = 1._wp
      heat_sum_eg ( jsb_hori ) = ( - 99._wp )
      chill_days_eg ( jsb_hori ) = 0._wp
     END IF
    END IF
   END IF
   IF ( abs ( growth_phase_eg ( jsb_hori ) ) < epsilon ( 1.0_wp ) ) THEN
    IF ( lai ( jsb_hori ) < phenoparam % all % laiseed ) THEN
     lai ( jsb_hori ) = min ( laimax_corrected ( jsb_hori ) - epsilon ( 1.0_wp&
      ) , phenoparam % all % laiseed )
    END IF
    IF ( lai ( jsb_hori ) > phenoparam % all % lai_negligible ) THEN
     grrate = phenoparam % eg % growthrate
     lai ( jsb_hori ) = get_letitgrow ( timestep_in_days , lai ( jsb_hori ) ,&
      grrate , 0.0_wp , laimax_corrected ( jsb_hori ) )
    END IF
   ELSE
    xtmp = phenoparam % eg % shedrate_vegetative
    lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
     xtmp )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_evergreen_phenology

 SUBROUTINE calc_crop_phenology ( day_of_year , is_newday , timestep_in_days ,&
  npp_pot_rate_ca , lat , pseudo_soil_temp , laimax_corrected ,&
  previous_day_npp_pot_rate_ca , w_soil_filling , previous_day_temp_min ,&
  previous_day_temp_max , lctlib_specificleafarea_c , heat_sum_crp ,&
  heat_sum_winter , growth_phase_crp , lai )
  USE mo_pheno_constants , ONLY: phenoparam

  INTEGER , INTENT(IN) :: day_of_year
  LOGICAL , INTENT(IN) :: is_newday
  REAL ( KIND= wp ) , INTENT(IN) :: timestep_in_days
  REAL ( KIND= wp ) , INTENT(IN) :: lctlib_specificleafarea_c
  REAL ( KIND= wp ) , INTENT(IN) :: npp_pot_rate_ca ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: lat ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: laimax_corrected ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: previous_day_npp_pot_rate_ca ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: pseudo_soil_temp ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: w_soil_filling ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: previous_day_temp_min ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: previous_day_temp_max ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: heat_sum_crp ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: heat_sum_winter ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: growth_phase_crp ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: lai ( : )
  REAL ( KIND= wp ) :: xtmp
  REAL ( KIND= wp ) :: grrate
  REAL ( KIND= wp ) :: temp_min_max
  LOGICAL :: npp_positive
  INTEGER :: jsb_hori

!$acc data &
!$acc present(npp_pot_rate_ca,lat,laimax_corrected &
!$acc ,previous_day_npp_pot_rate_ca,pseudo_soil_temp,w_soil_filling &
!$acc ,previous_day_temp_min,previous_day_temp_max,heat_sum_crp &
!$acc ,heat_sum_winter,growth_phase_crp,lai)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( previous_day_temp_max , 1 ) , 1
   IF ( previous_day_npp_pot_rate_ca ( jsb_hori ) > ( - 1.e-10_wp ) ) THEN
    npp_positive = .TRUE.
   ELSE
    npp_positive = .FALSE.
   END IF
   IF ( is_newday ) THEN
    IF ( lat ( jsb_hori ) >= 0._wp .AND. day_of_year == 70 .OR. lat ( jsb_hori&
     ) < 0._wp .AND. day_of_year == 252 ) THEN
     heat_sum_crp ( jsb_hori ) = 0._wp
     IF ( abs ( growth_phase_crp ( jsb_hori ) ) < 1.5_wp ) THEN
      growth_phase_crp ( jsb_hori ) = ( - 1._wp )
     END IF
    END IF
    IF ( lat ( jsb_hori ) >= 0._wp .AND. day_of_year == 172 .OR. lat (&
     jsb_hori ) < 0._wp .AND. day_of_year == 354 ) THEN
     IF ( abs ( growth_phase_crp ( jsb_hori ) ) < 1.5_wp .AND. heat_sum_crp (&
      jsb_hori ) > phenoparam % crp % heat_sum_harvest ) THEN
      heat_sum_crp ( jsb_hori ) = 0._wp
      growth_phase_crp ( jsb_hori ) = 1._wp
     ELSE
      IF ( growth_phase_crp ( jsb_hori ) > 1.5_wp ) THEN
       growth_phase_crp ( jsb_hori ) = ( - 2._wp )
      END IF
     END IF
     IF ( heat_sum_winter ( jsb_hori ) > 0._wp ) THEN
      heat_sum_winter ( jsb_hori ) = ( - 1._wp * heat_sum_winter ( jsb_hori )&
       )
     END IF
    END IF
    IF ( lat ( jsb_hori ) >= 0._wp .AND. day_of_year == 289 .OR. lat (&
     jsb_hori ) < 0._wp .AND. day_of_year == 105 ) THEN
     IF ( abs ( growth_phase_crp ( jsb_hori ) ) > 1.5_wp .AND. abs (&
      heat_sum_winter ( jsb_hori ) ) < phenoparam % crp % heat_sum_harvest&
      .AND. heat_sum_crp ( jsb_hori ) > phenoparam % crp % heat_sum_harvest )&
      THEN
      growth_phase_crp ( jsb_hori ) = ( - 1._wp )
     END IF
     IF ( abs ( growth_phase_crp ( jsb_hori ) ) > 1.5_wp .AND. heat_sum_crp (&
      jsb_hori ) > 2.0_wp * phenoparam % crp % heat_sum_harvest ) THEN
      growth_phase_crp ( jsb_hori ) = ( - 1._wp )
     END IF
     IF ( abs ( growth_phase_crp ( jsb_hori ) - 1._wp ) > 0.5_wp .AND.&
      heat_sum_crp ( jsb_hori ) < phenoparam % crp % heat_sum_harvest .AND.&
      abs ( heat_sum_winter ( jsb_hori ) ) > phenoparam % crp %&
      heat_sum_harvest .AND. pseudo_soil_temp ( jsb_hori ) > phenoparam % crp&
      % crit_temp ) THEN
      growth_phase_crp ( jsb_hori ) = ( - 2._wp )
     END IF
     heat_sum_winter ( jsb_hori ) = epsilon ( 1._wp )
    END IF
    temp_min_max = ( previous_day_temp_min ( jsb_hori ) +&
     previous_day_temp_max ( jsb_hori ) ) / 2.0_wp
    IF ( temp_min_max > phenoparam % crp % gdd_temp .AND. npp_positive ) THEN
     heat_sum_crp ( jsb_hori ) = heat_sum_crp ( jsb_hori ) + temp_min_max -&
      phenoparam % crp % gdd_temp
     IF ( heat_sum_winter ( jsb_hori ) > 0._wp ) THEN
      heat_sum_winter ( jsb_hori ) = heat_sum_winter ( jsb_hori ) +&
       temp_min_max - phenoparam % crp % gdd_temp
     END IF
    END IF
    IF ( pseudo_soil_temp ( jsb_hori ) > phenoparam % crp % crit_temp ) THEN
     IF ( abs ( growth_phase_crp ( jsb_hori ) ) < 1.5_wp ) THEN
      IF ( heat_sum_crp ( jsb_hori ) > phenoparam % crp % heat_sum_harvest )&
       THEN
       growth_phase_crp ( jsb_hori ) = ( - 1._wp )
      ELSE
       IF ( growth_phase_crp ( jsb_hori ) < 0.5_wp ) THEN
        growth_phase_crp ( jsb_hori ) = 0._wp
       END IF
      END IF
     ELSE
      IF ( heat_sum_winter ( jsb_hori ) > phenoparam % crp % heat_sum_harvest&
       .OR. heat_sum_winter ( jsb_hori ) < 0._wp ) THEN
       growth_phase_crp ( jsb_hori ) = ( - 2._wp )
      ELSE
       growth_phase_crp ( jsb_hori ) = 2._wp
      END IF
     END IF
    END IF
   END IF
   IF ( abs ( lat ( jsb_hori ) ) > 30._wp ) THEN
    IF ( growth_phase_crp ( jsb_hori ) > ( - 0.5_wp ) .AND. w_soil_filling (&
     jsb_hori ) > phenoparam % all % wilt_point + 1e-3_wp ) THEN
     IF ( pseudo_soil_temp ( jsb_hori ) > phenoparam % crp % crit_temp ) THEN
      IF ( lai ( jsb_hori ) < phenoparam % all % laiseed .AND. w_soil_filling&
       ( jsb_hori ) > phenoparam % crp % sproud ) THEN
       lai ( jsb_hori ) = min ( laimax_corrected ( jsb_hori ) - epsilon (&
        1.0_wp ) , phenoparam % all % laiseed )
      END IF
      IF ( previous_day_npp_pot_rate_ca ( jsb_hori ) > 0._wp ) THEN
       grrate = phenoparam % crp % leafalloc_fract * lctlib_specificleafarea_c&
        * max ( 0.0_wp , npp_pot_rate_ca ( jsb_hori ) ) * 86400.0_wp / lai (&
        jsb_hori )
       lai ( jsb_hori ) = get_letitgrow ( timestep_in_days , lai ( jsb_hori )&
        , grrate , 0.0_wp , laimax_corrected ( jsb_hori ) )
      ELSE
       xtmp = phenoparam % crp % shedrate_growth
       lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
        xtmp )
      END IF
     ELSE
      xtmp = phenoparam % crp % shedrate_growth
      lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
       xtmp )
     END IF
    ELSE
     xtmp = phenoparam % crp % shedrate_rest
     lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
      xtmp )
    END IF
   ELSE
    IF ( pseudo_soil_temp ( jsb_hori ) > phenoparam % crp % crit_temp .AND.&
     w_soil_filling ( jsb_hori ) > phenoparam % all % wilt_point + 1e-3_wp )&
     THEN
     IF ( lai ( jsb_hori ) < phenoparam % all % laiseed .AND. w_soil_filling (&
      jsb_hori ) > phenoparam % crp % sproud ) THEN
      lai ( jsb_hori ) = min ( laimax_corrected ( jsb_hori ) - epsilon (&
       1.0_wp ) , phenoparam % all % laiseed )
     END IF
     IF ( previous_day_npp_pot_rate_ca ( jsb_hori ) > 0._wp ) THEN
      grrate = phenoparam % crp % leafalloc_fract * lctlib_specificleafarea_c&
       * max ( 0.0_wp , npp_pot_rate_ca ( jsb_hori ) ) * 86400.0_wp / lai (&
       jsb_hori )
      lai ( jsb_hori ) = get_letitgrow ( timestep_in_days , lai ( jsb_hori ) ,&
       grrate , 0.0_wp , laimax_corrected ( jsb_hori ) )
     ELSE
      xtmp = phenoparam % crp % shedrate_growth
      lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
       xtmp )
     END IF
    ELSE
     xtmp = phenoparam % crp % shedrate_rest
     lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
      xtmp )
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_crop_phenology

 SUBROUTINE calc_raingreen_phenology ( timestep_in_days , w_soil_filling ,&
  laimax_corrected , previous_day_npp_pot_rate_ca , lai )
  USE mo_pheno_constants , ONLY: phenoparam

  REAL ( KIND= wp ) , INTENT(IN) :: timestep_in_days
  REAL ( KIND= wp ) , INTENT(IN) :: w_soil_filling ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: laimax_corrected ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: previous_day_npp_pot_rate_ca ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: lai ( : )
  REAL ( KIND= wp ) :: xtmp
  REAL ( KIND= wp ) :: grrate
  REAL ( KIND= wp ) :: average_filling
  INTEGER :: jsb_hori

!$acc data &
!$acc present(w_soil_filling,laimax_corrected,previous_day_npp_pot_rate_ca &
!$acc ,lai)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( laimax_corrected , 1 ) , 1
   average_filling = w_soil_filling ( jsb_hori )
   IF ( average_filling > phenoparam % all % wilt_point + 1e-3_wp ) THEN
    IF ( lai ( jsb_hori ) < phenoparam % all % laiseed .AND. average_filling >&
     phenoparam % rg % bucketfill_leafout ) THEN
     lai ( jsb_hori ) = min ( laimax_corrected ( jsb_hori ) - epsilon ( 1.0 )&
      , phenoparam % all % laiseed )
    END IF
    IF ( previous_day_npp_pot_rate_ca ( jsb_hori ) > 0._wp ) THEN
     grrate = phenoparam % rg % growthrate
     xtmp = get_shedrate_rg ( average_filling )
     lai ( jsb_hori ) = get_letitgrow ( timestep_in_days , lai ( jsb_hori ) ,&
      grrate , xtmp , laimax_corrected ( jsb_hori ) )
    END IF
   ELSE
    xtmp = phenoparam % rg % shedrate_dryseason
    lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
     xtmp )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_raingreen_phenology

 SUBROUTINE calc_grass_phenology ( t_air_in_celcius , timestep_in_days ,&
  w_soil_filling , laimax_corrected , previous_day_npp_pot_rate_ca , lai )
  USE mo_pheno_constants , ONLY: phenoparam

  REAL ( KIND= wp ) , INTENT(IN) :: timestep_in_days
  REAL ( KIND= wp ) , INTENT(IN) :: t_air_in_celcius ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: w_soil_filling ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: laimax_corrected ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: previous_day_npp_pot_rate_ca ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: lai ( : )
  REAL ( KIND= wp ) :: xtmp
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_air_in_celcius,w_soil_filling,laimax_corrected &
!$acc ,previous_day_npp_pot_rate_ca,lai)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( laimax_corrected , 1 ) , 1
   IF ( t_air_in_celcius ( jsb_hori ) > phenoparam % grs % crit_temp .AND.&
    w_soil_filling ( jsb_hori ) > phenoparam % all % wilt_point + 1e-3_wp )&
    THEN
    IF ( lai ( jsb_hori ) < phenoparam % all % laiseed ) THEN
     lai ( jsb_hori ) = min ( laimax_corrected ( jsb_hori ) - epsilon ( 1.0 )&
      , phenoparam % all % laiseed )
    END IF
    IF ( previous_day_npp_pot_rate_ca ( jsb_hori ) > 0._wp ) THEN
     xtmp = phenoparam % grs % growthrate
     lai ( jsb_hori ) = get_letitgrow ( timestep_in_days , lai ( jsb_hori ) ,&
      xtmp , 0.0_wp , laimax_corrected ( jsb_hori ) )
    ELSE
     xtmp = phenoparam % grs % shedrate_growth
     lai ( jsb_hori ) = get_letitgrow ( timestep_in_days , lai ( jsb_hori ) ,&
      0.0_wp , xtmp , laimax_corrected ( jsb_hori ) )
    END IF
   ELSE
    xtmp = phenoparam % grs % shedrate_dryseason
    lai ( jsb_hori ) = get_letitdie ( timestep_in_days , lai ( jsb_hori ) ,&
     xtmp )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_grass_phenology

 SUBROUTINE derive_maxlai_according_to_allometric_relationships (&
  lctlib_alpha_nr_ind , lctlib_beta_nr_ind , lctlib_alpha_leaf ,&
  lctlib_beta_leaf , lctlib_clumpinessfactor , lctlib_specificleafarea_c ,&
  c_woods , current_max_green , veg_carbon_at_max_green ,&
  number_of_individuals , biomass_per_individual , veg_fract_correction ,&
  maxlai_allom )
  USE mo_pheno_constants , ONLY: phenoparam
  USE mo_carbon_constants , ONLY: molarmassc_kg

  REAL ( KIND= wp ) , INTENT(IN) :: lctlib_alpha_nr_ind
  REAL ( KIND= wp ) , INTENT(IN) :: lctlib_beta_nr_ind
  REAL ( KIND= wp ) , INTENT(IN) :: lctlib_alpha_leaf
  REAL ( KIND= wp ) , INTENT(IN) :: lctlib_beta_leaf
  REAL ( KIND= wp ) , INTENT(IN) :: lctlib_clumpinessfactor
  REAL ( KIND= wp ) , INTENT(IN) :: lctlib_specificleafarea_c
  REAL ( KIND= wp ) , INTENT(IN) :: c_woods ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: current_max_green ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: veg_carbon_at_max_green ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: number_of_individuals ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: biomass_per_individual ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: veg_fract_correction ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: maxlai_allom ( : )
  REAL ( KIND= wp ) :: kgperha
  REAL ( KIND= wp ) :: leafc_per_individual
  REAL ( KIND= wp ) :: leafc
  INTEGER :: jsb_hori

!$acc data &
!$acc present(c_woods,current_max_green,veg_carbon_at_max_green &
!$acc ,number_of_individuals,biomass_per_individual,veg_fract_correction &
!$acc ,maxlai_allom)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( biomass_per_individual , 1 ) , 1
   kgperha = phenoparam % fr % kgctokg * molarmassc_kg * phenoparam % fr %&
    m2toha
   IF ( veg_carbon_at_max_green ( jsb_hori ) .gt. phenoparam % fr % min_c&
    .AND. c_woods ( jsb_hori ) .gt. phenoparam % fr % min_c ) THEN
    number_of_individuals ( jsb_hori ) = exp ( min ( ( log (&
     veg_carbon_at_max_green ( jsb_hori ) * kgperha ) - lctlib_alpha_nr_ind )&
     / lctlib_beta_nr_ind , phenoparam % fr % log_maxind ) )
    biomass_per_individual ( jsb_hori ) = max ( veg_carbon_at_max_green (&
     jsb_hori ) * kgperha / number_of_individuals ( jsb_hori ) , phenoparam %&
     fr % min_zero )
    leafc_per_individual = exp ( lctlib_alpha_leaf + lctlib_beta_leaf * log (&
     biomass_per_individual ( jsb_hori ) ) )
    leafc = leafc_per_individual * number_of_individuals ( jsb_hori ) /&
     kgperha
    maxlai_allom ( jsb_hori ) = max ( leafc * lctlib_specificleafarea_c ,&
     phenoparam % fr % min_maxlai )
   ELSE
    maxlai_allom ( jsb_hori ) = phenoparam % fr % min_maxlai
   END IF
   veg_fract_correction ( jsb_hori ) = 1.0_wp - exp ( ( - maxlai_allom (&
    jsb_hori ) / lctlib_clumpinessfactor ) )
   veg_carbon_at_max_green ( jsb_hori ) = 0.0_wp
   current_max_green ( jsb_hori ) = 0.0_wp
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE derive_maxlai_according_to_allometric_relationships

 SUBROUTINE track_max_green_pool ( c_green , c_woods , c_reserve ,&
  current_max_green , veg_carbon_at_max_green )

  REAL ( KIND= wp ) , INTENT(IN) :: c_green ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: c_woods ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: c_reserve ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: current_max_green ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: veg_carbon_at_max_green ( : )
  INTEGER :: jsb_hori

!$acc data &
!$acc present(c_green,c_woods,c_reserve,current_max_green &
!$acc ,veg_carbon_at_max_green)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( c_green , 1 ) , 1
   IF ( current_max_green ( jsb_hori ) + epsilon ( 1.0_wp ) .lt. c_green (&
    jsb_hori ) ) THEN
    current_max_green ( jsb_hori ) = c_green ( jsb_hori )
    veg_carbon_at_max_green ( jsb_hori ) = c_green ( jsb_hori ) + c_woods (&
     jsb_hori ) + c_reserve ( jsb_hori )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE track_max_green_pool

 FUNCTION get_letitdie ( timestep_in_days , x , p )
  REAL ( KIND= wp ) :: get_letitdie
  REAL ( KIND= wp ) , INTENT(IN) :: timestep_in_days
  REAL ( KIND= wp ) , INTENT(IN) :: x
  REAL ( KIND= wp ) , INTENT(IN) :: p

!$acc routine seq
  get_letitdie = exp ( ( - p * timestep_in_days ) ) * x
 END FUNCTION get_letitdie

 FUNCTION get_letitgrow ( timestep_in_days , x , k , p , z )
  REAL ( KIND= wp ) :: get_letitgrow
  REAL ( KIND= wp ) , INTENT(IN) :: timestep_in_days
  REAL ( KIND= wp ) , INTENT(IN) :: x
  REAL ( KIND= wp ) , INTENT(IN) :: k
  REAL ( KIND= wp ) , INTENT(IN) :: p
  REAL ( KIND= wp ) , INTENT(IN) :: z
  REAL ( KIND= wp ) :: hlp1
  REAL ( KIND= wp ) :: hlp2
  REAL ( KIND= wp ) :: numerator
  REAL ( KIND= wp ) :: denominator

!$acc routine seq
  hlp1 = k - p
  hlp2 = k * x / z
  numerator = hlp1
  denominator = hlp2 + exp ( ( - hlp1 * timestep_in_days ) ) * ( hlp1 - hlp2 )
  IF ( abs ( abs ( numerator ) - abs ( denominator ) ) <= 100._wp * epsilon (&
   1.0_wp ) ) THEN
   get_letitgrow = x
  ELSE
   get_letitgrow = x * numerator / denominator
  END IF
 END FUNCTION get_letitgrow

 FUNCTION get_shedrate_rg ( bucket_filling )
  USE mo_pheno_constants , ONLY: phenoparam
  REAL ( KIND= wp ) :: get_shedrate_rg
  REAL ( KIND= wp ) , INTENT(IN) :: bucket_filling
  REAL ( KIND= wp ) :: shedrate
  REAL ( KIND= wp ) :: slope

!$acc routine seq
  slope = phenoparam % rg % growthrate / ( phenoparam % rg %&
   bucketfill_critical - phenoparam % all % wilt_point )
  shedrate = phenoparam % rg % growthrate + phenoparam % rg % shedrate_aging -&
   slope * ( bucket_filling - phenoparam % all % wilt_point )
  get_shedrate_rg = min ( max ( phenoparam % rg % shedrate_aging , shedrate )&
   , phenoparam % rg % growthrate + phenoparam % rg % shedrate_aging )
 END FUNCTION get_shedrate_rg

 FUNCTION get_foliage_projected_cover ( fpc_max , lai , clumpiness )
  REAL ( KIND= wp ) :: get_foliage_projected_cover
  REAL ( KIND= wp ) , INTENT(IN) :: fpc_max
  REAL ( KIND= wp ) , INTENT(IN) :: lai
  REAL ( KIND= wp ) , INTENT(IN) :: clumpiness

!$acc routine seq
  get_foliage_projected_cover = fpc_max * ( 1._wp - exp ( ( - lai / clumpiness&
   ) ) )
 END FUNCTION get_foliage_projected_cover

END MODULE mo_pheno_process

