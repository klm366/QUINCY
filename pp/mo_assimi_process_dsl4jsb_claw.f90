MODULE mo_assimi_process
 USE mo_kind , ONLY: wp
 PUBLIC :: calc_assimilation_scaling_factors
 PUBLIC :: calc_assimilation_waterunlimited
 PUBLIC :: calc_assimilation_waterlimited
 PUBLIC :: calc_npp_pot_rate
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_assimi_process"

CONTAINS
 SUBROUTINE calc_assimilation_scaling_factors ( declination , canopy_bound_lai&
  , lai , lat , assimi_scaling_fact_cl )
  USE mo_rad_constants , ONLY: lailimit
  USE mo_math_constants , ONLY: pi

  REAL ( KIND= wp ) , INTENT(IN) :: declination
  REAL ( KIND= wp ) , INTENT(IN) :: canopy_bound_lai
  REAL ( KIND= wp ) , INTENT(IN) :: lai ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: lat ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: assimi_scaling_fact_cl ( : )
  REAL ( KIND= wp ) :: k12
  REAL ( KIND= wp ) :: cos_zenith_noon
  INTEGER :: jsb_hori

!$acc data present(lai,lat,assimi_scaling_fact_cl)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( lai , 1 ) , 1
   assimi_scaling_fact_cl ( jsb_hori ) = 1._wp
   cos_zenith_noon = cos ( declination ) * cos ( pi * lat ( jsb_hori ) /&
    180._wp ) + sin ( declination ) * sin ( pi * lat ( jsb_hori ) / 180._wp )
   k12 = 0.5_wp / max ( cos_zenith_noon , 1.e-3_wp )
   IF ( lai ( jsb_hori ) >= lailimit ) THEN
    assimi_scaling_fact_cl ( jsb_hori ) = max ( 1.e-10_wp , exp ( ( - k12 *&
     canopy_bound_lai * lai ( jsb_hori ) ) ) )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_assimilation_scaling_factors

 SUBROUTINE calc_assimilation_waterunlimited ( c4flag , carboxrate ,&
  etransport , t_air , press_srf , par_down_mol , apar_per_lai_cl ,&
  co2_mol_mixing_ratio , n_scaling_factors , canopy_cond_cl )
  USE mo_assimi_constants
  USE mo_jsb_physical_constants , ONLY: tmelt , r_gas => argas

  LOGICAL , INTENT(IN) :: c4flag
  REAL ( KIND= wp ) , INTENT(IN) :: carboxrate
  REAL ( KIND= wp ) , INTENT(IN) :: etransport
  REAL ( KIND= wp ) , INTENT(IN) :: t_air ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: press_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: par_down_mol ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: apar_per_lai_cl ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: co2_mol_mixing_ratio ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: n_scaling_factors ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: canopy_cond_cl ( : )
  REAL ( KIND= wp ) :: co2_conc_leaf_cl
  REAL ( KIND= wp ) :: t0
  REAL ( KIND= wp ) :: t1
  REAL ( KIND= wp ) :: kc
  REAL ( KIND= wp ) :: ko
  REAL ( KIND= wp ) :: gam
  REAL ( KIND= wp ) :: hitinhib
  REAL ( KIND= wp ) :: darkinhib
  REAL ( KIND= wp ) :: vcmax
  REAL ( KIND= wp ) :: jmax
  REAL ( KIND= wp ) :: k
  REAL ( KIND= wp ) :: tc
  REAL ( KIND= wp ) :: gass
  REAL ( KIND= wp ) :: dark_resp
  REAL ( KIND= wp ) :: je
  REAL ( KIND= wp ) :: jc
  REAL ( KIND= wp ) :: j0
  REAL ( KIND= wp ) :: j1
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_process" //&
  ":calc_assimilation_waterunlimited"
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_air,press_srf,par_down_mol,apar_per_lai_cl &
!$acc ,co2_mol_mixing_ratio,n_scaling_factors,canopy_cond_cl)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( t_air , 1 ) , 1
   IF ( c4flag ) THEN
    co2_conc_leaf_cl = fci1c4 * co2_mol_mixing_ratio ( jsb_hori )
   ELSE
    co2_conc_leaf_cl = fci1c3 * co2_mol_mixing_ratio ( jsb_hori )
   END IF
   t1 = 25._wp + tmelt
   tc = t_air ( jsb_hori ) - tmelt
   t0 = t_air ( jsb_hori ) - t1
   IF ( ( .NOT. c4flag ) ) THEN
    kc = kc0 * exp ( ec / r_gas * t0 * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori )&
     ** ( ( - 1 ) ) )
    ko = ko0 * exp ( eo / r_gas * t0 * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori )&
     ** ( ( - 1 ) ) )
    gam = max ( 1.7e-6_wp * tc , 0._wp )
    vcmax = carboxrate * n_scaling_factors ( jsb_hori ) * exp ( ev / r_gas *&
     t0 * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) )
    jmax = etransport * n_scaling_factors ( jsb_hori ) * tc / 25._wp
    jmax = max ( jmax , minofmaxcarboxrate )
    IF ( jmax .gt. minofmaxcarboxrate ) THEN
     j1 = alpha * apar_per_lai_cl ( jsb_hori ) * jmax / sqrt ( jmax ** ( 2 ) +&
      ( alpha * apar_per_lai_cl ( jsb_hori ) ) ** ( 2 ) )
    ELSE
     j1 = 0._wp
    END IF
    je = j1 * ( co2_conc_leaf_cl - gam ) / 4._wp / ( co2_conc_leaf_cl + 2._wp&
     * gam )
    jc = vcmax * ( co2_conc_leaf_cl - gam ) / ( co2_conc_leaf_cl + kc * (&
     1._wp + ox * ko ** ( ( - 1 ) ) ) )
    hitinhib = 1._wp / ( 1._wp + exp ( 1.3_wp * ( tc - 55._wp ) ) )
    gass = min ( je , jc ) * hitinhib
    darkinhib = 0.5_wp + 0.5_wp * exp ( ( - 2.e5_wp * max ( par_down_mol (&
     jsb_hori ) , 0.0_wp ) ) )
    dark_resp = frdc3 * carboxrate * n_scaling_factors ( jsb_hori ) * exp ( er&
     * t0 / r_gas * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) ) *&
     hitinhib * darkinhib
    canopy_cond_cl ( jsb_hori ) = max ( 1.6_wp * ( gass - dark_resp ) / (&
     co2_mol_mixing_ratio ( jsb_hori ) - co2_conc_leaf_cl ) * r_gas * t_air (&
     jsb_hori ) * press_srf ( jsb_hori ) ** ( ( - 1 ) ) , minstomaconductance&
     )
   ELSE
    IF ( c4flag ) THEN
     k = etransport * 1.e3_wp * n_scaling_factors ( jsb_hori ) * exp ( ek * t0&
      / r_gas * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) )
     vcmax = carboxrate * n_scaling_factors ( jsb_hori ) * exp ( ev / r_gas *&
      t0 * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) )
     darkinhib = 0.5_wp + 0.5_wp * exp ( ( - 2.e5_wp * max ( par_down_mol (&
      jsb_hori ) , 0.0_wp ) ) )
     hitinhib = 1._wp / ( 1._wp + exp ( 1.3_wp * ( tc - 55._wp ) ) )
     dark_resp = frdc4 * carboxrate * n_scaling_factors ( jsb_hori ) * exp (&
      er * t0 / r_gas * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) )&
      ) * hitinhib * darkinhib
     j0 = ( alc4 * apar_per_lai_cl ( jsb_hori ) + vcmax ) / 2._wp / theta
     je = j0 - sqrt ( j0 ** ( 2 ) - vcmax * alc4 * apar_per_lai_cl ( jsb_hori&
      ) / theta )
     jc = k * co2_conc_leaf_cl
     gass = min ( je , jc ) * hitinhib
     canopy_cond_cl ( jsb_hori ) = 1.6_wp * ( gass - dark_resp ) / (&
      co2_mol_mixing_ratio ( jsb_hori ) - co2_conc_leaf_cl ) * r_gas * t_air (&
      jsb_hori ) * press_srf ( jsb_hori ) ** ( ( - 1 ) )
     canopy_cond_cl ( jsb_hori ) = max ( canopy_cond_cl ( jsb_hori ) ,&
      minstomaconductance )
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_assimilation_waterunlimited

 SUBROUTINE calc_assimilation_waterlimited ( c4flag , carboxrate , etransport&
  , t_air , press_srf , par_down_mol , apar_per_lai_cl , co2_mol_mixing_ratio&
  , n_scaling_factors , co2_conc_leaf_cl , canopy_cond_cl ,&
  gross_assimilation_cl , dark_respiration_cl , carbox_rate_max_cl ,&
  e_transport_rate_max_cl , carbox_rate_cl , e_transport_rate_cl )
  USE mo_assimi_constants
  USE mo_physical_constants , ONLY: tmelt , r_gas => argas

  LOGICAL , INTENT(IN) :: c4flag
  REAL ( KIND= wp ) , INTENT(IN) :: carboxrate
  REAL ( KIND= wp ) , INTENT(IN) :: etransport
  REAL ( KIND= wp ) , INTENT(IN) :: t_air ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: press_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: par_down_mol ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: apar_per_lai_cl ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: co2_mol_mixing_ratio ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: n_scaling_factors ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: canopy_cond_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: co2_conc_leaf_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: gross_assimilation_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: dark_respiration_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: carbox_rate_max_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: e_transport_rate_max_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: carbox_rate_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: e_transport_rate_cl ( : )
  REAL ( KIND= wp ) :: t0
  REAL ( KIND= wp ) :: t1
  REAL ( KIND= wp ) :: kc
  REAL ( KIND= wp ) :: ko
  REAL ( KIND= wp ) :: gam
  REAL ( KIND= wp ) :: hitinhib
  REAL ( KIND= wp ) :: darkinhib
  REAL ( KIND= wp ) :: k1
  REAL ( KIND= wp ) :: k2
  REAL ( KIND= wp ) :: vcmax
  REAL ( KIND= wp ) :: jmax
  REAL ( KIND= wp ) :: k
  REAL ( KIND= wp ) :: tc
  REAL ( KIND= wp ) :: gass
  REAL ( KIND= wp ) :: dark_resp
  REAL ( KIND= wp ) :: je
  REAL ( KIND= wp ) :: jc
  REAL ( KIND= wp ) :: j0
  REAL ( KIND= wp ) :: j1
  REAL ( KIND= wp ) :: g0
  REAL ( KIND= wp ) :: w1
  REAL ( KIND= wp ) :: w2
  REAL ( KIND= wp ) :: b
  REAL ( KIND= wp ) :: c
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_process" //&
  ":calc_assimilation_waterlimited"
  INTEGER :: jsb_hori

!$acc data &
!$acc present(t_air,press_srf,par_down_mol,apar_per_lai_cl &
!$acc ,co2_mol_mixing_ratio,n_scaling_factors,canopy_cond_cl,co2_conc_leaf_cl &
!$acc ,gross_assimilation_cl,dark_respiration_cl,carbox_rate_max_cl &
!$acc ,e_transport_rate_max_cl,carbox_rate_cl,e_transport_rate_cl)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( apar_per_lai_cl , 1 ) , 1
   t1 = 25._wp + tmelt
   tc = t_air ( jsb_hori ) - tmelt
   t0 = t_air ( jsb_hori ) - t1
   IF ( ( .NOT. c4flag ) ) THEN
    kc = kc0 * exp ( ec / r_gas * t0 * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori )&
     ** ( ( - 1 ) ) )
    ko = ko0 * exp ( eo / r_gas * t0 * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori )&
     ** ( ( - 1 ) ) )
    gam = max ( 1.7e-6_wp * tc , 0._wp )
    vcmax = carboxrate * n_scaling_factors ( jsb_hori ) * exp ( ev * t0 /&
     r_gas * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) )
    darkinhib = 0.5_wp + 0.5_wp * exp ( ( - 2.e5_wp * max ( par_down_mol (&
     jsb_hori ) , 0.0_wp ) ) )
    hitinhib = 1._wp / ( 1._wp + exp ( 1.3_wp * ( tc - 55._wp ) ) )
    dark_resp = frdc3 * carboxrate * n_scaling_factors ( jsb_hori ) * exp ( er&
     * t0 / r_gas * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) ) *&
     hitinhib * darkinhib
    jmax = etransport * n_scaling_factors ( jsb_hori ) * tc / 25._wp
    jmax = max ( jmax , minofmaxcarboxrate )
    IF ( jmax .gt. minofmaxcarboxrate ) THEN
     j1 = alpha * apar_per_lai_cl ( jsb_hori ) * jmax / sqrt ( jmax ** ( 2 ) +&
      ( alpha * apar_per_lai_cl ( jsb_hori ) ) ** ( 2 ) )
    ELSE
     j1 = 0._wp
    END IF
    k1 = 2._wp * gam
    w1 = j1 / 4._wp
    w2 = vcmax
    k2 = kc * ( 1._wp + ox / ko )
    g0 = canopy_cond_cl ( jsb_hori ) / 1.6_wp / r_gas * t_air ( jsb_hori ) **&
     ( ( - 1 ) ) * press_srf ( jsb_hori )
    b = dark_resp + w1 + g0 * ( co2_mol_mixing_ratio ( jsb_hori ) + k1 )
    c = w1 * g0 * ( co2_mol_mixing_ratio ( jsb_hori ) - gam ) + w1 * dark_resp
    IF ( jmax .gt. minofmaxcarboxrate ) THEN
     je = b / 2._wp - sqrt ( max ( b ** ( 2 ) / 4._wp - c , 0._wp ) )
    ELSE
     je = 0._wp
    END IF
    b = dark_resp + w2 + g0 * ( co2_mol_mixing_ratio ( jsb_hori ) + k2 )
    c = w2 * g0 * ( co2_mol_mixing_ratio ( jsb_hori ) - gam ) + w2 * dark_resp
    jc = b / 2._wp - sqrt ( max ( b ** ( 2 ) / 4._wp - c , 0._wp ) )
    gass = min ( je , jc ) * hitinhib
    co2_conc_leaf_cl ( jsb_hori ) = co2_mol_mixing_ratio ( jsb_hori ) - max (&
     ( gass - dark_resp ) / max ( g0 , 1.e-6_wp ) , 0._wp )
   ELSE
    k = etransport * 1.e3_wp * n_scaling_factors ( jsb_hori ) * exp ( ek * t0&
     / r_gas * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) )
    vcmax = carboxrate * n_scaling_factors ( jsb_hori ) * exp ( ev * t0 /&
     r_gas * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) )
    darkinhib = 0.5_wp + 0.5_wp * exp ( ( - 2.e5_wp * max ( par_down_mol (&
     jsb_hori ) , 0.0_wp ) ) )
    hitinhib = 1._wp / ( 1._wp + exp ( 1.3_wp * ( tc - 55._wp ) ) )
    dark_resp = frdc4 * carboxrate * n_scaling_factors ( jsb_hori ) * exp ( er&
     * t0 / r_gas * t1 ** ( ( - 1 ) ) * t_air ( jsb_hori ) ** ( ( - 1 ) ) ) *&
     hitinhib * darkinhib
    g0 = canopy_cond_cl ( jsb_hori ) / 1.6_wp / r_gas * t_air ( jsb_hori ) **&
     ( ( - 1 ) ) * press_srf ( jsb_hori )
    j0 = ( alc4 * apar_per_lai_cl ( jsb_hori ) + vcmax ) / 2._wp / theta
    je = j0 - sqrt ( j0 ** ( 2 ) - vcmax * alc4 * apar_per_lai_cl ( jsb_hori )&
     / theta )
    jc = ( g0 * co2_mol_mixing_ratio ( jsb_hori ) + dark_resp ) / ( 1._wp + g0&
     / k )
    gass = min ( je , jc ) * hitinhib
    co2_conc_leaf_cl ( jsb_hori ) = co2_mol_mixing_ratio ( jsb_hori ) - max (&
     ( gass - dark_resp ) / max ( g0 , 1.e-6_wp ) , 0._wp )
    jmax = 0._wp
   END IF
   gross_assimilation_cl ( jsb_hori ) = gass
   dark_respiration_cl ( jsb_hori ) = dark_resp
   carbox_rate_max_cl ( jsb_hori ) = vcmax
   e_transport_rate_max_cl ( jsb_hori ) = jmax
   carbox_rate_cl ( jsb_hori ) = jc
   e_transport_rate_cl ( jsb_hori ) = je
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_assimilation_waterlimited

 SUBROUTINE calc_npp_pot_rate ( gross_assimilation_ca , dark_respiration_ca ,&
  npp_pot_rate_ca )
  USE mo_assimi_constants , ONLY: f_aut_leaf , ccost

  REAL ( KIND= wp ) , INTENT(IN) :: gross_assimilation_ca ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: dark_respiration_ca ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: npp_pot_rate_ca ( : )
  CHARACTER ( LEN= * ) , PARAMETER :: routine ="mo_assimi_process" //&
  ":calc_NPP_pot_rate"
  REAL ( KIND= wp ) , PARAMETER :: hlp = ( ccost - 1._wp ) / ccost
  REAL ( KIND= wp ) :: maintenancerespiration
  REAL ( KIND= wp ) :: growthrespiration
  INTEGER :: jsb_hori

!$acc data present(gross_assimilation_ca,dark_respiration_ca,npp_pot_rate_ca)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( dark_respiration_ca , 1 ) , 1
   maintenancerespiration = dark_respiration_ca ( jsb_hori ) / f_aut_leaf
   growthrespiration = max ( 0._wp , hlp * ( gross_assimilation_ca ( jsb_hori&
    ) - maintenancerespiration ) )
   npp_pot_rate_ca ( jsb_hori ) = gross_assimilation_ca ( jsb_hori ) -&
    maintenancerespiration - growthrespiration
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_npp_pot_rate

END MODULE mo_assimi_process

