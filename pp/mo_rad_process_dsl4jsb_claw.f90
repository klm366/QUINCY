MODULE mo_rad_process
 USE mo_kind , ONLY: wp
 PRIVATE :: has_minimal_radiation_2
 PRIVATE :: has_minimal_radiation_1
 PUBLIC :: calc_radiation_surface_net
 PUBLIC :: calc_par
 PUBLIC :: calc_alb_lwtr
 PUBLIC :: calc_alb_lice
 PUBLIC :: calc_glacier_albedo
 PUBLIC :: has_minimal_radiation
 PUBLIC :: calc_soil_albedo
 PUBLIC :: calc_snow_albedo
 PUBLIC :: calc_sky_view_fractions
 PUBLIC :: merge_albedos_of_vegtile
 PUBLIC :: get_surface_albedo_simple
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_rad_process"
 INTERFACE has_minimal_radiation
  MODULE PROCEDURE has_minimal_radiation_1
  MODULE PROCEDURE has_minimal_radiation_2
 END INTERFACE

CONTAINS
 SUBROUTINE calc_radiation_surface_net ( swvis_down , swnir_down , alb_vis ,&
  alb_nir , lw_down , t , rad_net , swvis_net , swnir_net , sw_net , lw_net )
  USE mo_phy_schemes , ONLY: lwnet_from_lwdown

  REAL ( KIND= wp ) , INTENT(IN) :: swvis_down ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: swnir_down ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_vis ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_nir ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: lw_down ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: rad_net ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(OUT) :: swvis_net ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(OUT) :: swnir_net ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(OUT) :: sw_net ( : )
  REAL ( KIND= wp ) , OPTIONAL , INTENT(OUT) :: lw_net ( : )
  REAL ( KIND= wp ) :: zswvis_net
  REAL ( KIND= wp ) :: zswnir_net
  REAL ( KIND= wp ) :: zsw_net
  REAL ( KIND= wp ) :: zlw_net
  INTEGER :: jsb_hori

!$acc data &
!$acc present(swvis_down,swnir_down,alb_vis,alb_nir,lw_down,t,rad_net &
!$acc ,swvis_net,swnir_net,sw_net,lw_net)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( swnir_down , 1 ) , 1
   zswvis_net = swvis_down ( jsb_hori ) * ( 1._wp - alb_vis ( jsb_hori ) )
   zswnir_net = swnir_down ( jsb_hori ) * ( 1._wp - alb_nir ( jsb_hori ) )
   zsw_net = zswvis_net + zswnir_net
   zlw_net = lwnet_from_lwdown ( lw_down ( jsb_hori ) , t ( jsb_hori ) )
   rad_net ( jsb_hori ) = zsw_net + zlw_net
   IF ( present ( swvis_net ) ) THEN
    swvis_net ( jsb_hori ) = zswvis_net
   END IF
   IF ( present ( swnir_net ) ) THEN
    swnir_net ( jsb_hori ) = zswnir_net
   END IF
   IF ( present ( sw_net ) ) THEN
    sw_net ( jsb_hori ) = zsw_net
   END IF
   IF ( present ( lw_net ) ) THEN
    lw_net ( jsb_hori ) = zlw_net
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_radiation_surface_net

 SUBROUTINE calc_par ( ncanopy , icanopy , use_alb_veg_simple , par ,&
  alb_vis_soil , lai , cos_zenith_angle , fract_par_direct , canopy_bound_lai&
  , canopy_bound_lai_delta , b4_layer_above , par_down_mol ,&
  soil_reflectivity_par , fapar_cl , lai_cl , apar_per_lai_cl )
  USE mo_rad_constants , ONLY: epar , soilreflectivityparmin , fcmax , fcmin ,&
   zenithminpar , laimin , lailimit

  INTEGER , INTENT(IN) :: ncanopy
  INTEGER , INTENT(IN) :: icanopy
  LOGICAL , INTENT(IN) :: use_alb_veg_simple
  REAL ( KIND= wp ) , INTENT(IN) :: canopy_bound_lai
  REAL ( KIND= wp ) , INTENT(IN) :: canopy_bound_lai_delta
  REAL ( KIND= wp ) , INTENT(IN) :: par ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_vis_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: lai ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: cos_zenith_angle ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: fract_par_direct ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: b4_layer_above ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: par_down_mol ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: soil_reflectivity_par ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: fapar_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: lai_cl ( : )
  REAL ( KIND= wp ) , INTENT(OUT) :: apar_per_lai_cl ( : )
  REAL ( KIND= wp ) :: zh
  REAL ( KIND= wp ) :: k0
  REAL ( KIND= wp ) :: zp1
  REAL ( KIND= wp ) :: q0
  REAL ( KIND= wp ) :: q1
  REAL ( KIND= wp ) :: b0
  REAL ( KIND= wp ) :: b1
  REAL ( KIND= wp ) :: ekl
  REAL ( KIND= wp ) :: ehl
  REAL ( KIND= wp ) :: b4
  REAL ( KIND= wp ) :: fc
  REAL ( KIND= wp ) :: fshd
  REAL ( KIND= wp ) :: zp0
  REAL ( KIND= wp ) :: ekl0
  REAL ( KIND= wp ) :: ehl0
  REAL ( KIND= wp ) :: x0
  REAL ( KIND= wp ) :: x1
  REAL ( KIND= wp ) :: x2
  REAL ( KIND= wp ) :: f
  REAL ( KIND= wp ) , PARAMETER :: omega = 0.12_wp
  INTEGER :: jsb_hori

!$acc data &
!$acc present(par,alb_vis_soil,lai,cos_zenith_angle,fract_par_direct &
!$acc ,b4_layer_above,par_down_mol,soil_reflectivity_par,fapar_cl,lai_cl &
!$acc ,apar_per_lai_cl)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( par , 1 ) , 1
   par_down_mol ( jsb_hori ) = par ( jsb_hori ) / epar
   IF ( ( .NOT. use_alb_veg_simple ) ) THEN
    soil_reflectivity_par ( jsb_hori ) = alb_vis_soil ( jsb_hori )
   ELSE
    soil_reflectivity_par ( jsb_hori ) = max ( 0.92_wp * alb_vis_soil (&
     jsb_hori ) - 0.015_wp , soilreflectivityparmin )
   END IF
   fapar_cl ( jsb_hori ) = 0._wp
   apar_per_lai_cl ( jsb_hori ) = 0._wp
   lai_cl ( jsb_hori ) = lai ( jsb_hori ) * canopy_bound_lai_delta
   lai_cl ( jsb_hori ) = max ( lai_cl ( jsb_hori ) , laimin )
   IF ( lai ( jsb_hori ) < lailimit ) THEN
    fc = lai ( jsb_hori ) / lailimit * fcmax
   ELSE
    fc = fcmax
   END IF
   fc = max ( fc , fcmin )
   IF ( cos_zenith_angle ( jsb_hori ) .ge. zenithminpar ) THEN
    zh = sqrt ( 1._wp - omega )
    zp1 = ( 1._wp - omega / 2._wp + zh ) / omega * 2._wp
    zp0 = 1._wp / zp1
    k0 = 0.5_wp / cos_zenith_angle ( jsb_hori )
    IF ( k0 == zh ) THEN
     k0 = k0 + 1.e-12_wp
    END IF
    IF ( k0 == ( - zh ) ) THEN
     k0 = k0 + 1.e-12_wp
    END IF
    x0 = 1._wp - 4._wp * cos_zenith_angle ( jsb_hori ) ** ( 2 ) * zh ** ( 2 )
    q1 = ( 1._wp + 2._wp * cos_zenith_angle ( jsb_hori ) ) * omega / 2._wp /&
     x0
    q0 = ( 1._wp - 2._wp * cos_zenith_angle ( jsb_hori ) ) * omega / 2._wp /&
     x0
    fshd = max ( fc , fcmin )
    ekl = exp ( ( - k0 / fshd * lai ( jsb_hori ) ) )
    ehl = exp ( ( - zh / fshd * lai ( jsb_hori ) ) )
    x1 = ( zp1 - soil_reflectivity_par ( jsb_hori ) ) / ehl
    x0 = ( zp0 - soil_reflectivity_par ( jsb_hori ) ) * ehl
    x2 = ( soil_reflectivity_par ( jsb_hori ) * ( 1._wp - q1 ) - q0 ) * ekl
    f = 1._wp - fract_par_direct ( jsb_hori ) + q1 * fract_par_direct (&
     jsb_hori )
    b1 = ( x2 * fract_par_direct ( jsb_hori ) - f * x0 ) / ( x1 - x0 )
    b0 = x2 * fract_par_direct ( jsb_hori ) / ( x0 - x1 ) + f / ( 1.0_wp - x0&
     / x1 )
    b4 = ( 1._wp - q0 - q1 ) * fract_par_direct ( jsb_hori ) + ( 1._wp - zp1 )&
     * b1 + ( 1._wp - zp0 ) * b0
   END IF
   IF ( b4_layer_above ( jsb_hori ) < 7777776._wp ) THEN
    b4 = b4_layer_above ( jsb_hori )
   END IF
   IF ( cos_zenith_angle ( jsb_hori ) .ge. zenithminpar ) THEN
    IF ( icanopy < ncanopy ) THEN
     zp0 = 1._wp / zp1
     ekl0 = exp ( ( - k0 / fshd * canopy_bound_lai * lai ( jsb_hori ) ) )
     ehl0 = exp ( ( - zh / fshd * canopy_bound_lai * lai ( jsb_hori ) ) )
     x0 = ( 1._wp - q0 - q1 ) * ekl0 * fract_par_direct ( jsb_hori ) + ( 1._wp&
      - zp1 ) * b1 / ehl0 + ( 1._wp - zp0 ) * b0 * ehl0
     fapar_cl ( jsb_hori ) = b4 - x0
     b4_layer_above ( jsb_hori ) = x0
    ELSE
     x0 = ( 1._wp - q0 - q1 ) * ekl * fract_par_direct ( jsb_hori ) + ( 1._wp&
      - zp1 ) * b1 / ehl + ( 1._wp - zp0 ) * b0 * ehl
     fapar_cl ( jsb_hori ) = b4_layer_above ( jsb_hori ) - x0
    END IF
   END IF
   IF ( cos_zenith_angle ( jsb_hori ) .ge. zenithminpar ) THEN
    fapar_cl ( jsb_hori ) = fapar_cl ( jsb_hori ) * fshd
   END IF
   apar_per_lai_cl ( jsb_hori ) = par_down_mol ( jsb_hori ) * fapar_cl (&
    jsb_hori ) / max ( lai_cl ( jsb_hori ) , 1.e-10_wp )
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_par

 SUBROUTINE calc_alb_lwtr ( l_day , alb_lwtr )
  USE mo_rad_constants , ONLY: albedolakewater

  LOGICAL , INTENT(IN) :: l_day ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_lwtr ( : )
  INTEGER :: jsb_hori

!$acc data present(l_day,alb_lwtr)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( l_day , 1 ) , 1
   IF ( ( .NOT. l_day ( jsb_hori ) ) ) THEN
!cdir CLAW: RETURN statement transformed for parallel region
   ELSE
    alb_lwtr ( jsb_hori ) = albedolakewater
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_alb_lwtr

 SUBROUTINE calc_alb_lice ( l_day , t , w_snow_lice , alb_lice )
  USE mo_rad_constants , ONLY: albedolakeicemin , albedolakeicemax ,&
   albedolakesnowmin , albedolakesnowmax
  USE mo_jsb_physical_constants , ONLY: tmelt

  LOGICAL , INTENT(IN) :: l_day ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: w_snow_lice ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_lice ( : )
  REAL ( KIND= wp ) :: t_upper_limit
  REAL ( KIND= wp ) :: alb_min
  REAL ( KIND= wp ) :: alb_max
  INTEGER :: jsb_hori

!$acc data present(l_day,t,w_snow_lice,alb_lice)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( w_snow_lice , 1 ) , 1
   IF ( ( .NOT. l_day ( jsb_hori ) ) ) THEN
!cdir CLAW: RETURN statement transformed for parallel region
   ELSE
    t_upper_limit = tmelt - 1._wp
    IF ( w_snow_lice ( jsb_hori ) > 0.01_wp ) THEN
     alb_min = albedolakesnowmin
     alb_max = albedolakesnowmax
    ELSE
     alb_min = albedolakeicemin
     alb_max = albedolakeicemax
    END IF
    IF ( t ( jsb_hori ) >= tmelt ) THEN
     alb_lice ( jsb_hori ) = alb_min
    ELSE
     IF ( t ( jsb_hori ) < t_upper_limit ) THEN
      alb_lice ( jsb_hori ) = alb_max
     ELSE
      alb_lice ( jsb_hori ) = alb_min + ( alb_max - alb_min ) / ( tmelt -&
       t_upper_limit ) * ( tmelt - t ( jsb_hori ) )
     END IF
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_alb_lice

 SUBROUTINE calc_glacier_albedo ( l_day , t , alb_vis , alb_nir )
  USE mo_rad_constants , ONLY: albedoglaciervismin , albedoglaciervismax ,&
   tempalbedoglaciermax , albedoglaciernirmin , albedoglaciernirmax
  USE mo_jsb_physical_constants , ONLY: tmelt

  LOGICAL , INTENT(IN) :: l_day ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_vis ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_nir ( : )
  INTEGER :: jsb_hori

!$acc data present(l_day,t,alb_vis,alb_nir)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( alb_vis , 1 ) , 1
   IF ( ( .NOT. l_day ( jsb_hori ) ) ) THEN
!cdir CLAW: RETURN statement transformed for parallel region
   ELSE
    IF ( t ( jsb_hori ) >= tmelt ) THEN
     alb_vis ( jsb_hori ) = 0.76
     alb_nir ( jsb_hori ) = 0.34
    ELSE
     IF ( t ( jsb_hori ) < tmelt - 5.0 ) THEN
      alb_vis ( jsb_hori ) = 0.89
      alb_nir ( jsb_hori ) = 0.72
     ELSE
      alb_vis ( jsb_hori ) = 0.76 + ( tmelt - t ( jsb_hori ) ) * ( 0.89 - 0.76&
       ) / 5.0
      alb_nir ( jsb_hori ) = 0.34 + ( tmelt - t ( jsb_hori ) ) * ( 0.72 - 0.34&
       ) / 5.0
     END IF
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_glacier_albedo

 ELEMENTAL SUBROUTINE calc_soil_albedo ( l_day , use_alb_mineralsoil_const ,&
  use_alb_soil_organic_int , use_alb_soil_litter , c_ag_sum_1 , c_bg_sum ,&
  specificleafarea_c , albedolittervis , alb_background , alb_vis_mineralsoil&
  , alb_nir_mineralsoil , alb_vis_soil , alb_nir_soil )
  USE mo_rad_constants , ONLY: albedosoilconstvis , albedosoilconstnir

  LOGICAL , INTENT(IN) :: l_day
  LOGICAL , INTENT(IN) :: use_alb_mineralsoil_const
  INTEGER , INTENT(IN) :: use_alb_soil_organic_int
  LOGICAL , INTENT(IN) :: use_alb_soil_litter
  REAL ( KIND= wp ) , INTENT(IN) :: c_ag_sum_1
  REAL ( KIND= wp ) , INTENT(IN) :: c_bg_sum
  REAL ( KIND= wp ) , INTENT(IN) :: specificleafarea_c
  REAL ( KIND= wp ) , INTENT(IN) :: albedolittervis
  REAL ( KIND= wp ) , INTENT(IN) :: alb_background
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_vis_mineralsoil
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_nir_mineralsoil
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_vis_soil
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_nir_soil
  REAL ( KIND= wp ) , PARAMETER :: scalfact_soc_to_alb_vis = 3.e-4
  REAL ( KIND= wp ) , PARAMETER :: scalfact_soc_to_alb_nir = 3.e-4
  REAL ( KIND= wp ) , PARAMETER :: soil_c_min = 500.
  REAL ( KIND= wp ) :: alb_vis_soil_nolitter
  REAL ( KIND= wp ) :: alb_nir_soil_nolitter
  REAL ( KIND= wp ) :: litter_area
  REAL ( KIND= wp ) :: litter_view_fract
  REAL ( KIND= wp ) :: alb_vis_litter
  REAL ( KIND= wp ) :: alb_nir_litter

  IF ( ( .NOT. l_day ) ) THEN
   RETURN
  END IF
  IF ( use_alb_mineralsoil_const ) THEN
   alb_vis_mineralsoil = 0.2
   alb_nir_mineralsoil = 0.2
  ELSE
   alb_vis_mineralsoil = alb_background
   alb_nir_mineralsoil = alb_background
  END IF
  IF ( use_alb_soil_organic_int .eq. 1 ) THEN
   alb_vis_soil_nolitter = alb_vis_mineralsoil - scalfact_soc_to_alb_vis * min&
    ( c_bg_sum , soil_c_min )
   alb_nir_soil_nolitter = alb_nir_mineralsoil - scalfact_soc_to_alb_nir * min&
    ( c_bg_sum , soil_c_min )
  ELSE
   IF ( use_alb_soil_organic_int .eq. 2 ) THEN
    alb_vis_soil_nolitter = max ( 6._wp * alb_vis_mineralsoil / ( 6._wp + log&
     ( 1._wp + c_bg_sum ) ) , 0.1_wp )
    alb_nir_soil_nolitter = max ( 6._wp * alb_nir_mineralsoil / ( 6._wp + log&
     ( 1._wp + c_bg_sum ) ) , 0.2_wp )
   ELSE
    alb_vis_soil_nolitter = alb_vis_mineralsoil
    alb_nir_soil_nolitter = alb_nir_mineralsoil
   END IF
  END IF
  IF ( use_alb_soil_litter ) THEN
   litter_area = c_ag_sum_1 * specificleafarea_c
   alb_vis_litter = litter_area * albedolittervis
   IF ( litter_area <= epsilon ( 1._wp ) ) THEN
    alb_vis_litter = 0.1_wp
    alb_nir_litter = 0.2_wp
   END IF
   litter_view_fract = 1.0_wp - exp ( ( - 0.5_wp * litter_area ) )
   alb_vis_soil = ( 1.0_wp - litter_view_fract ) * alb_vis_soil_nolitter +&
    litter_view_fract * alb_vis_litter
   alb_nir_soil = ( 1.0_wp - litter_view_fract ) * alb_nir_soil_nolitter +&
    litter_view_fract * alb_nir_litter
  ELSE
   alb_vis_soil = alb_vis_soil_nolitter
   alb_nir_soil = alb_nir_soil_nolitter
  END IF
 END SUBROUTINE calc_soil_albedo

 SUBROUTINE calc_snow_albedo ( l_day , fract_snow_soil , snow_age ,&
  cos_zenith_angle , albedo_age_weight , t_srf , alb_vis_snow , alb_nir_snow )
  USE mo_rad_constants , ONLY: zenithanglefactor , albedosnowvismax_age ,&
   albedosnownirmax_age , albedosnowvisage , albedosnownirage ,&
   albedosnowangle , tempalbedosnowmax , albedosnowvismin_temp ,&
   albedosnownirmin_temp , albedosnowvismax_temp , albedosnownirmax_temp
  USE mo_jsb_physical_constants , ONLY: tmelt

  LOGICAL , INTENT(IN) :: l_day ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: fract_snow_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: snow_age ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: cos_zenith_angle ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: albedo_age_weight
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_vis_snow ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_nir_snow ( : )
  REAL ( KIND= wp ) :: alb_vis_snow_age
  REAL ( KIND= wp ) :: alb_nir_snow_age
  REAL ( KIND= wp ) :: alb_vis_snow_temp
  REAL ( KIND= wp ) :: alb_nir_snow_temp
  REAL ( KIND= wp ) :: snow_age_factor
  REAL ( KIND= wp ) :: snow_albedo_angle_factor
  INTEGER :: jsb_hori

!$acc data &
!$acc present(l_day,fract_snow_soil,snow_age,cos_zenith_angle,t_srf &
!$acc ,alb_vis_snow,alb_nir_snow)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( cos_zenith_angle , 1 ) , 1
   alb_vis_snow_age = 0.0_wp
   alb_nir_snow_age = 0.0_wp
   alb_vis_snow_temp = 0.0_wp
   alb_nir_snow_temp = 0.0_wp
   alb_vis_snow ( jsb_hori ) = 0.0_wp
   alb_nir_snow ( jsb_hori ) = 0.0_wp
   IF ( ( .NOT. l_day ( jsb_hori ) ) ) THEN
!cdir CLAW: RETURN statement transformed for parallel region
   ELSE
    IF ( fract_snow_soil ( jsb_hori ) > epsilon ( 1._wp ) ) THEN
     snow_age_factor = snow_age ( jsb_hori ) / ( 1._wp + snow_age ( jsb_hori )&
      )
     IF ( cos_zenith_angle ( jsb_hori ) < 0.5_wp ) THEN
      snow_albedo_angle_factor = ( ( 1.0_wp + 2. ) / ( 1.0_wp + 2.0_wp * 2. *&
       cos_zenith_angle ( jsb_hori ) ) - 1.0_wp ) / 2.
     ELSE
      snow_albedo_angle_factor = 0._wp
     END IF
     alb_vis_snow_age = 0.90 * ( 1._wp - 0.15 * snow_age_factor )
     alb_vis_snow_age = alb_vis_snow_age + 0.4 * snow_albedo_angle_factor * (&
      1._wp - alb_vis_snow_age )
     alb_nir_snow_age = 0.60 * ( 1._wp - 0.5 * snow_age_factor )
     alb_nir_snow_age = alb_nir_snow_age + 0.4 * snow_albedo_angle_factor * (&
      1._wp - alb_nir_snow_age )
    ELSE
     alb_vis_snow_age = 0._wp
     alb_nir_snow_age = 0._wp
    END IF
    IF ( t_srf ( jsb_hori ) >= tmelt ) THEN
     alb_vis_snow_temp = 0.52
     alb_nir_snow_temp = 0.3
    ELSE
     IF ( t_srf ( jsb_hori ) < tmelt - 5.0 ) THEN
      alb_vis_snow_temp = 0.90
      alb_nir_snow_temp = 0.65
     ELSE
      alb_vis_snow_temp = 0.52 + ( tmelt - t_srf ( jsb_hori ) ) * ( 0.90 -&
       0.52 ) / 5.0
      alb_nir_snow_temp = 0.3 + ( tmelt - t_srf ( jsb_hori ) ) * ( 0.65 - 0.3&
       ) / 5.0
     END IF
    END IF
    alb_vis_snow ( jsb_hori ) = albedo_age_weight * alb_vis_snow_age + ( 1._wp&
     - albedo_age_weight ) * alb_vis_snow_temp
    alb_nir_snow ( jsb_hori ) = albedo_age_weight * alb_nir_snow_age + ( 1._wp&
     - albedo_age_weight ) * alb_nir_snow_temp
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_snow_albedo

 SUBROUTINE calc_sky_view_fractions ( l_day , lai , fract_fpc_max ,&
  cos_zenith_angle , stemarea , sky_view_fract , sky_view_fract_stem )
  USE mo_rad_constants , ONLY: skyviewfactor

  LOGICAL , INTENT(IN) :: l_day ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: lai ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: fract_fpc_max ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: cos_zenith_angle ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: stemarea
  REAL ( KIND= wp ) , INTENT(INOUT) :: sky_view_fract ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: sky_view_fract_stem ( : )
  INTEGER :: jsb_hori

!$acc data &
!$acc present(l_day,lai,fract_fpc_max,cos_zenith_angle,sky_view_fract &
!$acc ,sky_view_fract_stem)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( cos_zenith_angle , 1 ) , 1
   IF ( ( .NOT. l_day ( jsb_hori ) ) ) THEN
!cdir CLAW: RETURN statement transformed for parallel region
   ELSE
    sky_view_fract ( jsb_hori ) = 1._wp - fract_fpc_max ( jsb_hori ) * ( 1._wp&
     - exp ( ( - 0.5 * lai ( jsb_hori ) ) ) )
    sky_view_fract_stem ( jsb_hori ) = sky_view_fract ( jsb_hori ) - ( 1._wp -&
     fract_fpc_max ( jsb_hori ) * ( 1._wp - exp ( ( - 0.5 * ( lai ( jsb_hori )&
     + stemarea ) * ( 2.0_wp - cos_zenith_angle ( jsb_hori ) ) ) ) ) )
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE calc_sky_view_fractions

 SUBROUTINE merge_albedos_of_vegtile ( l_day , l_forest , sky_view_fract ,&
  sky_view_fract_stem , albedocanopysnow , fract_snow_soil , fract_snow_can ,&
  alb_vis_snow , alb_nir_snow , alb_vis_soil , alb_nir_soil , alb_vis_can ,&
  alb_nir_can , alb_vis , alb_nir )

  LOGICAL , INTENT(IN) :: l_day ( : )
  LOGICAL , INTENT(IN) :: l_forest
  REAL ( KIND= wp ) , INTENT(IN) :: sky_view_fract ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: sky_view_fract_stem ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: albedocanopysnow
  REAL ( KIND= wp ) , INTENT(IN) :: fract_snow_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: fract_snow_can ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_vis_snow ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_nir_snow ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_vis_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_nir_soil ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_vis_can ( : )
  REAL ( KIND= wp ) , INTENT(IN) :: alb_nir_can ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_vis ( : )
  REAL ( KIND= wp ) , INTENT(INOUT) :: alb_nir ( : )
  REAL ( KIND= wp ) :: alb_vis_nostems
  REAL ( KIND= wp ) :: alb_nir_nostems
  INTEGER :: jsb_hori

!$acc data &
!$acc present(l_day,sky_view_fract,sky_view_fract_stem,fract_snow_soil &
!$acc ,fract_snow_can,alb_vis_snow,alb_nir_snow,alb_vis_soil,alb_nir_soil &
!$acc ,alb_vis_can,alb_nir_can,alb_vis,alb_nir)
!$acc parallel
!$acc loop gang vector
  DO jsb_hori = 1 , size ( alb_vis , 1 ) , 1
   IF ( ( .NOT. l_day ( jsb_hori ) ) ) THEN
!cdir CLAW: RETURN statement transformed for parallel region
   ELSE
    alb_vis_nostems = ( 1.0_wp - sky_view_fract ( jsb_hori ) ) * alb_vis_can (&
     jsb_hori ) + sky_view_fract ( jsb_hori ) * alb_vis_soil ( jsb_hori )
    alb_nir_nostems = ( 1.0_wp - sky_view_fract ( jsb_hori ) ) * alb_nir_can (&
     jsb_hori ) + sky_view_fract ( jsb_hori ) * alb_nir_soil ( jsb_hori )
    IF ( l_forest ) THEN
     alb_vis ( jsb_hori ) = max ( ( sky_view_fract ( jsb_hori ) -&
      sky_view_fract_stem ( jsb_hori ) ) * ( fract_snow_soil ( jsb_hori ) *&
      alb_vis_snow ( jsb_hori ) + ( 1.0_wp - fract_snow_soil ( jsb_hori ) ) *&
      alb_vis_soil ( jsb_hori ) ) + sky_view_fract_stem ( jsb_hori ) * (&
      fract_snow_can ( jsb_hori ) * albedocanopysnow + ( 1.0_wp -&
      fract_snow_can ( jsb_hori ) ) * alb_vis_soil ( jsb_hori ) ) + ( 1.0_wp -&
      sky_view_fract ( jsb_hori ) ) * ( fract_snow_can ( jsb_hori ) *&
      albedocanopysnow + ( 1.0_wp - fract_snow_can ( jsb_hori ) ) *&
      alb_vis_can ( jsb_hori ) ) , alb_vis_nostems )
     alb_nir ( jsb_hori ) = max ( ( sky_view_fract ( jsb_hori ) -&
      sky_view_fract_stem ( jsb_hori ) ) * ( fract_snow_soil ( jsb_hori ) *&
      alb_nir_snow ( jsb_hori ) + ( 1.0_wp - fract_snow_soil ( jsb_hori ) ) *&
      alb_nir_soil ( jsb_hori ) ) + sky_view_fract_stem ( jsb_hori ) * (&
      fract_snow_can ( jsb_hori ) * albedocanopysnow + ( 1.0_wp -&
      fract_snow_can ( jsb_hori ) ) * alb_nir_soil ( jsb_hori ) ) + ( 1.0_wp -&
      sky_view_fract ( jsb_hori ) ) * ( fract_snow_can ( jsb_hori ) *&
      albedocanopysnow + ( 1.0_wp - fract_snow_can ( jsb_hori ) ) *&
      alb_nir_can ( jsb_hori ) ) , alb_nir_nostems )
    ELSE
     alb_vis ( jsb_hori ) = max ( fract_snow_soil ( jsb_hori ) * alb_vis_snow&
      ( jsb_hori ) + ( 1.0_wp - fract_snow_soil ( jsb_hori ) ) *&
      alb_vis_nostems , alb_vis_nostems )
     alb_nir ( jsb_hori ) = max ( fract_snow_soil ( jsb_hori ) * alb_nir_snow&
      ( jsb_hori ) + ( 1.0_wp - fract_snow_soil ( jsb_hori ) ) *&
      alb_nir_nostems , alb_nir_nostems )
    END IF
   END IF
  END DO
!$acc end parallel
!$acc end data
 END SUBROUTINE merge_albedos_of_vegtile

 FUNCTION get_surface_albedo_simple ( t_srf , fract_forest , grnd_snow_fract ,&
  fract_snow_can , lai , alb_background )
  USE mo_jsb_physical_constants , ONLY: tmelt
  USE mo_rad_constants , ONLY: albedocanopysnow_simple , albedosnowmin_simple&
   , albedosnowmax_simple
  REAL ( KIND= wp ) :: get_surface_albedo_simple
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf
  REAL ( KIND= wp ) , INTENT(IN) :: fract_forest
  REAL ( KIND= wp ) , INTENT(IN) :: grnd_snow_fract
  REAL ( KIND= wp ) , INTENT(IN) :: fract_snow_can
  REAL ( KIND= wp ) , INTENT(IN) :: lai
  REAL ( KIND= wp ) , INTENT(IN) :: alb_background
  REAL ( KIND= wp ) :: sky_view_fract
  REAL ( KIND= wp ) :: forest_view_fract
  REAL ( KIND= wp ) :: albedo_snow_grnd
  REAL ( KIND= wp ) :: min_temp_snow_albedo
  REAL ( KIND= wp ) :: albedo_grnd
  REAL ( KIND= wp ) :: albedo_can

!$acc routine seq
  min_temp_snow_albedo = tmelt - 5._wp
  IF ( t_srf >= tmelt ) THEN
   albedo_snow_grnd = 0.4
  ELSE
   IF ( t_srf < min_temp_snow_albedo ) THEN
    albedo_snow_grnd = 0.8
   ELSE
    albedo_snow_grnd = 0.4 + ( tmelt - t_srf ) * ( 0.8 - 0.4 ) / ( tmelt -&
     min_temp_snow_albedo )
   END IF
  END IF
  albedo_grnd = grnd_snow_fract * albedo_snow_grnd + ( 1._wp - grnd_snow_fract&
   ) * alb_background
  albedo_can = fract_snow_can * 0.2 + ( 1._wp - fract_snow_can ) *&
   alb_background
  sky_view_fract = exp ( ( - 1.0_wp * max ( lai , 2._wp ) ) )
  forest_view_fract = fract_forest * ( 1._wp - sky_view_fract )
  get_surface_albedo_simple = max ( alb_background , ( 1._wp -&
   forest_view_fract ) * albedo_grnd + forest_view_fract * albedo_can )
 END FUNCTION get_surface_albedo_simple

 FUNCTION has_minimal_radiation_2 ( swvis , swnir ) RESULT(l_day)
  USE mo_rad_constants , ONLY: minradiation
  REAL ( KIND= wp ) , INTENT(IN) :: swvis
  REAL ( KIND= wp ) , INTENT(IN) :: swnir
  LOGICAL :: l_day

!$acc routine seq
  l_day = swvis + swnir > 1.0e-09
 END FUNCTION has_minimal_radiation_2

 FUNCTION has_minimal_radiation_1 ( sw ) RESULT(l_day)
  USE mo_rad_constants , ONLY: minradiation
  REAL ( KIND= wp ) , INTENT(IN) :: sw
  LOGICAL :: l_day

!$acc routine seq
  l_day = sw > 1.0e-09
 END FUNCTION has_minimal_radiation_1

END MODULE mo_rad_process

