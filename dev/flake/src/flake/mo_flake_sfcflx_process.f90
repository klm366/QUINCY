!>
!! @brief Contains the routines for the processes to compute the surface fluxes from waterbodies
!!
!! @author
!!  Dmitrii Mironov, DWD
!!  Zoé Rehder, MPI-M
!!
!! @par Revision History
!! First version by Dmitrii Mironov (2005/11/17)
!! Adapted for JSBACH by Zoé Rehder (2023/02/09)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!!  Description:
!!  The SfcFlx routine 
!!  where fluxes of momentum and of sensible and latent heat 
!!  at the air-water or air-ice (air-snow) interface are computed. 

MODULE mo_flake_SfcFlx_process
#ifndef __NO_JSBACH__
  USE mo_kind,              ONLY: wp
      
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                           &
    &   calc_momsenlat              , &
    &   calc_roughness

  CHARACTER(len=*), PARAMETER :: modname = 'mo_flake_SfcFlx_process'

  CONTAINS

  PURE ELEMENTAL SUBROUTINE calc_momsenlat (&
    &  height_u   , & ! in
    &  height_tq  , & ! in
    &  fetch      , & ! in
    &  U_a        , & ! in
    &  T_a        , & ! in
    &  q_a        , & ! in
    &  T_s        , & ! in
    &  P_a        , & ! in
    &  h_ice      , & ! in
    &  Q_momentum , & ! out
    &  Q_sensible , & ! out
    &  Q_latent   , & ! out
    &  Q_watvap   , & ! out
    &  z0u_sf     , & ! out
    &  z0t_sf     , & ! out
    &  u_star_st    & ! out
    )

    USE mo_atmland_util,              ONLY: calc_spec_humidity_sat
    USE mo_jsb_physical_constants,    ONLY: alf                     , & !< Latent heat of fusion [J kg^{-1}]
                                            alv                     , & !< Specific heat of evaporation [J kg^{-1}]
                                            rd                      , & !< Gas constant for dry air               [J/(K kg)]
                                            rv                      , & !< Gas constant for water vapor           [J/(K*kg)]
                                            cpd                     , & !< Specific heat of air at constant pressure     [J kg^{-1} K^{-1}]
                                            grav                        !< Acceleration due to gravity [m s^{-2}]

    USE mo_flake_SfcFlx_constants,  ONLY : kappa_q_a               , & !< Molecular diffusivity of air for water vapour [m^{2} s^{-1}]  
                                            c_free_conv             , & !< Constant in the expressions for fluxes in free convection 
                                            c_Karman                , & !< The von Karman constant 
                                            c_MO_u_exp              , & !< Constant of the MO theory (wind, exponent)
                                            c_MO_t_stab             , & !< Constant of the MO theory (wind, stable stratification)
                                            c_MO_q_stab             , & !< Constant of the MO theory (humidity, stable stratification)
                                            c_MO_u_stab             , & !< Constant of the MO theory (temperature, stable stratification)
                                            c_MO_t_conv             , & !< Constant of the MO theory (temperature, convection)
                                            c_MO_u_conv             , & !< Constant of the MO theory (wind, convection)                                     
                                            c_MO_q_conv             , & !< Constant of the MO theory (humidity, convection)
                                            c_MO_t_exp              , & !< Constant of the MO theory (temperature, exponent)
                                            c_MO_q_exp              , & !< Constant of the MO theory (humidity, exponent)
                                            Pr_neutral              , & !< Turbulent Prandtl number at neutral static stability 
                                            Sc_neutral                  !< Turbulent Schmidt number at neutral static stability

    USE mo_flake_parameters,        ONLY : c_accur_sf              , & !< A small number (accuracy)
                                            u_star_min_sf           , & !< Minimum value of friction velocity [m s^{-1}]
                                            u_wind_min_sf           , & !< Minimum wind speed [m s^{-1}]
                                            c_small                 , & !< A small number
                                            h_Ice_min                   !< Minimum ice thickness [m]

    USE mo_flake_constants,         ONLY : con_m                   , & !< Kinematic viscosity of dry air         [m^2/s]  
                                            con_h                       !< Scalar conductivity of dry air         [m^2/s]  

    !$claw model-data
    REAL(wp), INTENT(in) ::  height_u   , & ! Height where wind is measured [m]
                          &  height_tq  , & ! Height where temperature and humidity are measured [m]
                          &  fetch      , & ! Typical wind fetch [m]
                          &  U_a        , & ! Wind speed [m s^{-1}]
                          &  T_a        , & ! Air temperature [K]
                          &  q_a        , & ! Air specific humidity [-]
                          &  T_s        , & ! Surface temperature (water, ice or snow) [K]
                          &  P_a        , & ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]
                          &  h_ice           ! Ice thickness [m]

    REAL(wp),  INTENT(out) :: Q_momentum    , & ! Momentum flux [N m^{-2}]  
                            & Q_sensible    , & ! Sensible heat flux [W m^{-2}]  
                            & Q_latent      , & ! Laten heat flux [W m^{-2}]
                            & Q_watvap      , & ! Flux of water vapour [kg m^{-2} s^{-1}]
                            & z0u_sf        , & ! Roughness length with respect to wind velocity [m]
                            & z0t_sf        , & ! Roughness length with respect to potential temperature [m]
                            & u_star_st         ! Friction velocity [m s^{-1}]
    !$claw end model-data

    INTEGER, PARAMETER ::  n_iter_max     = 24  ! Maximum number of iterations 
    
    !  Local variables of type LOGICAL
    LOGICAL ::  l_conv_visc     , & ! Switch, TRUE = viscous free convection, the Nu=C Ra^(1/3) law is used
                l_conv_cbl          ! Switch, TRUE = CBL scale convective structures define surface fluxes 


    !  Local variables of type INTEGER
    INTEGER  :: i        , & ! Loop index
                n_iter       ! Number of iterations performed 

    !  Local variables of type REAL
    REAL (KIND = wp) :: rho_a                  , & ! Air density [kg m^{-3}]  
                        q_s                        ! Saturation specific humidity at T=T_s [-]

    !  Roughness lengths
    REAL (KIND = wp) :: z0q_sf                     ! Roughness length with respect to specific humidity [m]

    !  Local variables of type REAL
    REAL (KIND = wp) :: Q_mom_tur              , & ! Turbulent momentum flux [N m^{-2}]
                        Q_sen_tur              , & ! Turbulent sensible heat flux [W m^{-2}]  
                        Q_lat_tur              , & ! Turbulent laten heat flux [W m^{-2}]
                        Q_mom_mol              , & ! Molecular momentum flux [N m^{-2}]
                        Q_sen_mol              , & ! Molecular sensible heat flux [W m^{-2}]  
                        Q_lat_mol              , & ! Molecular laten heat flux [W m^{-2}]
                        Q_mom_con              , & ! Momentum flux in free convection [N m^{-2}]
                        Q_sen_con              , & ! Sensible heat flux in free convection [W m^{-2}]  
                        Q_lat_con                  ! Laten heat flux in free convection [W m^{-2}]

    !  Local variables of type REAL
    REAL (KIND = wp) :: par_conv_visc          , & ! Viscous convection stability parameter
                        par_conv_cbl           , & ! CBL convection stability parameter
                        c_z0u_fetch            , & ! Fetch-dependent Charnock parameter
                        U_a_thresh             , & ! Threshld value of the wind speed [m s^{-1}] 
                        u_star_thresh          , & ! Threshld value of friction velocity [m s^{-1}]
                        u_star_previter        , & ! Friction velocity from previous iteration [m s^{-1}]
                        u_star_n               , & ! Friction velocity at neutral stratification [m s^{-1}]
                        ! u_star_st              , & ! Friction velocity with due regard for stratification [m s^{-1}]
                        ZoL                    , & ! The z/L ratio, z=height_u
                        Ri                     , & ! Gradient Richardson number 
                        Ri_cr                  , & ! Critical value of Ri 
                        R_z                    , & ! Ratio of "height_tq" to "height_u"
                        Fun                    , & ! A function of generic variable "x"
                        Fun_prime              , & ! Derivative of "Fun" with respect to "x"
                        Delta                  , & ! Relative error 
                        psi_u                  , & ! The MO stability function for wind profile
                        psi_t                  , & ! The MO stability function for temperature profile
                        psi_q                      ! The MO stability function for specific humidity profile         
                        
    !  Derived thermodynamic parameters 
    REAL (KIND = wp), PARAMETER ::  Rd_o_Rv  = rd/rv                   , & ! Ratio of gas constants (Rd/Rv)
                                    alpha_q  = (1._wp-Rd_o_Rv)/Rd_o_Rv     ! Dimensionless ratio                            
    ! ---------------------------
    ! Go      
    ! All fluxes are positive when directed upwards.  
                        
    !  Compute saturation specific humidity and the air density at T=T_s    
    q_s       = calc_spec_humidity_sat(T_s, P_a)         ! Saturation specific humidity at T=T_s
    rho_a     = P_a/rd/T_s/(1._wp+(1._wp/Rd_o_Rv-1._wp)*q_s)     ! Air density at T_s and q_s (surface values)   
    
    !  Compute molecular fluxes of momentum and of sensible and latent heat
    !_dm The fluxes are in kinematic units
    Q_mom_mol = -con_m*U_a/height_u 
    Q_sen_mol = -con_h*(T_a-T_s)/height_tq    
    Q_lat_mol = -kappa_q_a*(q_a-q_s)/height_tq

    !  Compute fluxes in free convection
    par_conv_visc = (T_s-T_a)/T_s*SQRT(con_h) + (q_s-q_a)*alpha_q*SQRT(kappa_q_a)
    IF(par_conv_visc.GT.0._wp) THEN   ! Viscous convection takes place
      l_conv_visc = .TRUE.
      par_conv_visc = (par_conv_visc*grav/con_m)**(1._wp/3._wp)
      Q_sen_con = c_free_conv*SQRT(con_h)*par_conv_visc  
      Q_sen_con = Q_sen_con*(T_s-T_a)
      Q_lat_con = c_free_conv*SQRT(kappa_q_a)*par_conv_visc
      Q_lat_con = Q_lat_con*(q_s-q_a)
    ELSE                                  ! No viscous convection, set fluxes to zero
      l_conv_visc = .FALSE.
      Q_sen_con = 0._wp 
      Q_lat_con = 0._wp
    END IF
    Q_mom_con = 0._wp                 ! Momentum flux in free (viscous or CBL-scale) convection is zero  

    !  Compute turbulent fluxes

    R_z   = height_tq/height_u                        ! Ratio of "height_tq" to "height_u"
    Ri_cr = c_MO_t_stab/c_MO_u_stab**2*R_z            ! Critical Ri
    Ri    = grav*((T_a-T_s)/T_s+alpha_q*(q_a-q_s))/MAX(U_a,u_wind_min_sf)**2_wp
    Ri    = Ri*height_u/Pr_neutral                    ! Gradient Richardson number

    IF(U_a.LT.u_wind_min_sf.OR.Ri.GT.Ri_cr-c_small) THEN  ! Low wind or Ri>Ri_cr

      u_star_st = 0._wp                               ! Set turbulent fluxes to zero 
      Q_mom_tur = 0._wp                       
      Q_sen_tur = 0._wp   
      Q_lat_tur = 0._wp  

    ELSE                                              ! Compute turbulent fluxes using MO similarity

      ! Compute z/L, where z=height_u
      IF(Ri.GE.0._wp) THEN                            ! Stable stratification
        ZoL = SQRT(1._wp-4._wp*(c_MO_u_stab-R_z*c_MO_t_stab)*Ri)
        ZoL = ZoL - 1._wp + 2._wp*c_MO_u_stab*Ri
        ZoL = ZoL/2._wp/c_MO_u_stab/c_MO_u_stab/(Ri_cr-Ri)
      ELSE                                            ! Convection
        n_iter = 0_wp
        Delta = 1._wp                                 ! Set initial error to a large value (as compared to the accuracy)
        u_star_previter = Ri*MAX(1._wp, SQRT(R_z*c_MO_t_conv/c_MO_u_conv)) ! Initial guess for ZoL
        DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
          Fun = u_star_previter**2_wp*(c_MO_u_conv*u_star_previter-1._wp)  &
              + Ri**2_wp*(1._wp-R_z*c_MO_t_conv*u_star_previter)
          Fun_prime = 3._wp*c_MO_u_conv*u_star_previter**2_wp              &
                    - 2._wp*u_star_previter - R_z*c_MO_t_conv*Ri**2_wp
          ZoL = u_star_previter - Fun/Fun_prime
          Delta = ABS(ZoL-u_star_previter)/MAX(c_accur_sf, ABS(ZoL+u_star_previter))
          u_star_previter = ZoL
          n_iter = n_iter + 1_wp
        END DO 
      END IF
      !  Compute fetch-dependent Charnock parameter, use "u_star_min_sf"
      CALL calc_roughness (fetch, U_a, u_star_min_sf, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)

      !  Threshold value of wind speed 
      u_star_st = u_star_thresh
      CALL calc_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)

      IF(ZoL.GT.0._wp) THEN   ! MO function in stable stratification 
        psi_u = c_MO_u_stab*ZoL*(1._wp-MIN(z0u_sf/height_u, 1._wp))
      ELSE                        ! MO function in convection
        psi_t = (1._wp-c_MO_u_conv*ZoL)**c_MO_u_exp
        psi_q = (1._wp-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1._wp))**c_MO_u_exp
        psi_u = 2._wp*(ATAN(psi_t)-ATAN(psi_q))                  &
              + 2._wp*LOG((1._wp+psi_q)/(1._wp+psi_t))   &
              + LOG((1._wp+psi_q*psi_q)/(1._wp+psi_t*psi_t))   
      END IF 
      U_a_thresh = u_star_thresh/c_Karman*(LOG(height_u/z0u_sf)+psi_u)

      !  Compute friction velocity 
      n_iter = 0
      Delta = 1._wp                ! Set initial error to a large value (as compared to the accuracy)
      u_star_previter = u_star_thresh  ! Initial guess for friction velocity  
      IF(U_a.LE.U_a_thresh) THEN  ! Smooth surface
        DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
          CALL calc_roughness (fetch, U_a, MIN(u_star_thresh, u_star_previter), h_ice,   &
                                c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
          IF(ZoL.GE.0._wp) THEN  ! Stable stratification
            psi_u = c_MO_u_stab*ZoL*(1._wp-MIN(z0u_sf/height_u, 1._wp))
            Fun = LOG(height_u/z0u_sf) + psi_u
            Fun_prime = (Fun + 1._wp + c_MO_u_stab*ZoL*MIN(z0u_sf/height_u, 1._wp))/c_Karman
            Fun = Fun*u_star_previter/c_Karman - U_a
          ELSE                       ! Convection 
            psi_t = (1._wp-c_MO_u_conv*ZoL)**c_MO_u_exp
            psi_q = (1._wp-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1._wp))**c_MO_u_exp
            psi_u = 2._wp*(ATAN(psi_t)-ATAN(psi_q))                  &
                  + 2._wp*LOG((1._wp+psi_q)/(1._wp+psi_t))   &
                  + LOG((1._wp+psi_q*psi_q)/(1._wp+psi_t*psi_t))   
            Fun = LOG(height_u/z0u_sf) + psi_u
            Fun_prime = (Fun + 1._wp/psi_q)/c_Karman
            Fun = Fun*u_star_previter/c_Karman - U_a
          END IF
          u_star_st = u_star_previter - Fun/Fun_prime
          Delta = ABS((u_star_st-u_star_previter)/(u_star_st+u_star_previter))
          u_star_previter = u_star_st
          n_iter = n_iter + 1
        END DO 
      ELSE                        ! Rough surface 
        DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
          CALL calc_roughness (fetch, U_a, MAX(u_star_thresh, u_star_previter), h_ice,   &
                                c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
          IF(ZoL.GE.0._wp) THEN  ! Stable stratification
            psi_u = c_MO_u_stab*ZoL*(1._wp-MIN(z0u_sf/height_u, 1._wp))
            Fun = LOG(height_u/z0u_sf) + psi_u
            Fun_prime = (Fun - 2._wp - 2._wp*c_MO_u_stab*ZoL*MIN(z0u_sf/height_u, 1._wp))/c_Karman
            Fun = Fun*u_star_previter/c_Karman - U_a
          ELSE                       ! Convection 
            psi_t = (1._wp-c_MO_u_conv*ZoL)**c_MO_u_exp
            psi_q = (1._wp-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1._wp))**c_MO_u_exp
            psi_u = 2._wp*(ATAN(psi_t)-ATAN(psi_q))                  &
                  + 2._wp*LOG((1._wp+psi_q)/(1._wp+psi_t))   &
                  + LOG((1._wp+psi_q*psi_q)/(1._wp+psi_t*psi_t))   
            Fun = LOG(height_u/z0u_sf) + psi_u
            Fun_prime = (Fun - 2._wp/psi_q)/c_Karman
            Fun = Fun*u_star_previter/c_Karman - U_a
          END IF
          IF(h_ice.GE.h_Ice_min) THEN   ! No iteration is required for rough flow over ice
            u_star_st = c_Karman*U_a/MAX(c_small, LOG(height_u/z0u_sf)+psi_u)
            u_star_previter = u_star_st
          ELSE                              ! Iterate in case of open water
            u_star_st = u_star_previter - Fun/Fun_prime
          END IF
          Delta = ABS((u_star_st-u_star_previter)/(u_star_st+u_star_previter))
          u_star_previter = u_star_st
          n_iter = n_iter + 1
        END DO 
      END IF
      
      !  Momentum flux
      Q_mom_tur = -u_star_st*u_star_st

      !  Temperature and specific humidity fluxes
      CALL calc_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
      IF(ZoL.GE.0._wp) THEN   ! Stable stratification 
        psi_t = c_MO_t_stab*R_z*ZoL*(1._wp-MIN(z0t_sf/height_tq, 1._wp))
        psi_q = c_MO_q_stab*R_z*ZoL*(1._wp-MIN(z0q_sf/height_tq, 1._wp))

      ELSE                        ! Convection 
        psi_u = (1._wp-c_MO_t_conv*R_z*ZoL)**c_MO_t_exp
        psi_t = (1._wp-c_MO_t_conv*R_z*ZoL*MIN(z0t_sf/height_tq, 1._wp))**c_MO_t_exp
        psi_t = 2._wp*LOG((1._wp+psi_t)/(1._wp+psi_u))
        psi_u = (1._wp-c_MO_q_conv*R_z*ZoL)**c_MO_q_exp
        psi_q = (1._wp-c_MO_q_conv*R_z*ZoL*MIN(z0q_sf/height_tq, 1._wp))**c_MO_q_exp
        psi_q = 2._wp*LOG((1._wp+psi_q)/(1._wp+psi_u))

      END IF 
      Q_sen_tur = -(T_a-T_s)*u_star_st*c_Karman/Pr_neutral  &
                / MAX(c_small, LOG(height_tq/z0t_sf)+psi_t)
      Q_lat_tur = -(q_a-q_s)*u_star_st*c_Karman/Sc_neutral  &
                / MAX(c_small, LOG(height_tq/z0q_sf)+psi_q)
      
    END IF   ! Turbulent fluxes
  
    !  Decide between turbulent, molecular, and convective fluxes
    Q_momentum = MIN(Q_mom_tur, Q_mom_mol, Q_mom_con)  ! Momentum flux is negative          
    IF(l_conv_visc) THEN    ! Convection, take fluxes that are maximal in magnitude 
      IF(ABS(Q_sen_tur).GE.ABS(Q_sen_con)) THEN
        Q_sensible = Q_sen_tur
      ELSE
        Q_sensible = Q_sen_con
      END IF
      IF(ABS(Q_sensible).LT.ABS(Q_sen_mol)) THEN
        Q_sensible = Q_sen_mol
      END IF
      IF(ABS(Q_lat_tur).GE.ABS(Q_lat_con)) THEN
        Q_latent = Q_lat_tur
      ELSE
        Q_latent = Q_lat_con
      END IF
      IF(ABS(Q_latent).LT.ABS(Q_lat_mol)) THEN
        Q_latent = Q_lat_mol
      END IF
    ELSE                    ! Stable or neutral stratification, chose fluxes that are maximal in magnitude 
      IF(ABS(Q_sen_tur).GE.ABS(Q_sen_mol)) THEN 
        Q_sensible = Q_sen_tur
      ELSE 
        Q_sensible = Q_sen_mol    
      END IF
      IF(ABS(Q_lat_tur).GE.ABS(Q_lat_mol)) THEN 
        Q_latent = Q_lat_tur
      ELSE 
        Q_latent = Q_lat_mol  
      END IF
    END IF 

    !  Set output (notice that fluxes are no longer in kinematic units)

    Q_momentum = Q_momentum*rho_a 
    Q_sensible = Q_sensible*rho_a*cpd
    Q_watvap   = Q_latent*rho_a
    Q_latent = alv
    IF(h_ice.GE.h_Ice_min) Q_latent = Q_latent + alf   ! Add latent heat of fusion over ice
    Q_latent = Q_watvap*Q_latent

  END SUBROUTINE calc_momsenlat

  PURE ELEMENTAL SUBROUTINE calc_roughness(&
    &  fetch          , & ! in
    &  U_a            , & ! in
    &  u_star         , & ! in
    &  h_ice          , & ! in
    &  c_z0u_fetch    , & ! out
    &  u_star_thresh  , & ! out
    &  z0u            , & ! out
    &  z0t            , & ! out
    &  z0q              & ! out
    )

    ! Description:
    !
    !  Computes the water-surface or the ice-surface roughness lengths
    !  with respect to wind velocity, potential temperature and specific humidity.
    !
    !  The water-surface roughness lengths with respect to wind velocity is computed
    !  from the Charnock formula when the surface is aerodynamically rough.
    !  A simple empirical formulation is used to account for the dependence 
    !  of the Charnock parameter on the wind fetch. 
    !  When the flow is aerodynamically smooth, the roughness length with respect to 
    !  wind velocity is proportional to the depth of the viscous sub-layer.
    !  The water-surface roughness lengths for scalars are computed using the power-law 
    !  formulations in terms of the roughness Reynolds number (Zilitinkevich et al. 2001).
    !  The ice-surface aerodynamic roughness is taken to be constant.
    !  The ice-surface roughness lengths for scalars 
    !  are computed through the power-law formulations 
    !  in terms of the roughness Reynolds number (Andreas 2002).  
    
    USE mo_jsb_physical_constants,    ONLY: grav                        !< Acceleration due to gravity [m s^{-2}]
    USE mo_flake_constants,          ONLY: con_m                       !< Kinematic viscosity of dry air         [m^2/s]  
    USE mo_flake_parameters,         ONLY: c_accur_sf              , & !< A small number (accuracy)
                                            h_Ice_min                ,& !< Minimum ice thickness [m]
                                            u_wind_min_sf               !< Minimum wind speed [m s^{-1}]


    USE mo_flake_SfcFlx_constants,  ONLY : c_Karman                , & !< The von Karman constant 
                                            c_z0u_smooth            , & !< Constant in the expression for z0u (smooth flow) 
                                            c_z0u_rough             , & !< The Charnock constant in the expression for z0u (rough flow)
                                            c_z0u_rough_L           , & !< An increased Charnock constant (used as the upper limit)
                                            c_z0u_ftch_f            , & !< Factor in the expression for fetch-dependent Charnock parameter
                                            c_z0u_ftch_ex           , & !< Exponent in the expression for fetch-dependent Charnock parameter
                                            c_z0t_rough_1           , & !< Constant in the expression for z0t (factor) 
                                            c_z0t_rough_2           , & !< Constant in the expression for z0t (factor)
                                            c_z0t_rough_3           , & !< Constant in the expression for z0t (exponent) 
                                            c_z0q_rough_1           , & !< Constant in the expression for z0q (factor)
                                            c_z0q_rough_2           , & !< Constant in the expression for z0q (factor)
                                            c_z0q_rough_3           , & !< Constant in the expression for z0q (exponent)
                                            c_z0t_ice_b0s           , & !< Constant in the expression for z0t over ice
                                            c_z0t_ice_b0t           , & !< Constant in the expression for z0t over ice
                                            c_z0t_ice_b1t           , & !< Constant in the expression for z0t over ice
                                            c_z0t_ice_b0r           , & !< Constant in the expression for z0t over ice
                                            c_z0t_ice_b1r           , & !< Constant in the expression for z0t over ice
                                            c_z0t_ice_b2r           , & !< Constant in the expression for z0t over ice
                                            c_z0q_ice_b0s           , & !< Constant in the expression for z0q over ice
                                            c_z0q_ice_b0t           , & !< Constant in the expression for z0q over ice
                                            c_z0q_ice_b1t           , & !< Constant in the expression for z0q over ice
                                            c_z0q_ice_b0r           , & !< Constant in the expression for z0q over ice
                                            c_z0q_ice_b1r           , & !< Constant in the expression for z0q over ice
                                            c_z0q_ice_b2r           , & !< Constant in the expression for z0q over ice
                                            Re_z0s_ice_t            , & !< Threshold value of the surface Reynolds number 
                                            Pr_neutral              , & !< Turbulent Prandtl number at neutral static stability 
                                            Sc_neutral              , & !< Turbulent Schmidt number at neutral static stability
                                            z0u_ice_rough               !< Aerodynamic roughness of the ice surface [m] (rough flow)

    !$claw model-data
    REAL(wp), INTENT(in) :: fetch    , & !< Typical wind fetch [m]
                            U_a      , & !< Wind speed [m s^{-1}]
                            u_star   , & !< Friction velocity in the surface air layer [m s^{-1}]
                            h_ice        !< Ice thickness [m]
    
    REAL(wp), INTENT(out):: c_z0u_fetch    , & !< Fetch-dependent Charnock parameter
                            u_star_thresh  , & !< Threshold value of friction velocity [m s^{-1}]
                            z0u            , & !< Roughness length with respect to wind velocity [m]
                            z0t            , & !< Roughness length with respect to potential temperature [m]
                            z0q                !< Roughness length with respect to specific humidity [m]
    !$claw end model-data

    !  Local variables
    REAL (KIND = wp) :: Re_s          , & !< Surface Reynolds number 
                        Re_s_thresh       !< Threshold value of Re_s

    ! ---------------------------
    ! Go   
    IF(h_ice.LT.h_Ice_min) THEN  ! Water surface  

      ! The Charnock parameter as dependent on dimensionless fetch
      c_z0u_fetch = MAX(U_a, u_wind_min_sf)**2_wp/grav/fetch  ! Inverse dimensionless fetch
      c_z0u_fetch = c_z0u_rough + c_z0u_ftch_f*c_z0u_fetch**c_z0u_ftch_ex
      c_z0u_fetch = MIN(c_z0u_fetch, c_z0u_rough_L)                      ! Limit Charnock parameter
      
      ! Threshold value of friction velocity
      u_star_thresh = (c_z0u_smooth/c_z0u_fetch*grav*con_m)**(1._wp/3._wp)
      
      ! Surface Reynolds number and its threshold value
      Re_s = u_star**3_wp/con_m/grav
      Re_s_thresh = c_z0u_smooth/c_z0u_fetch
      
      ! Aerodynamic roughness
      IF(Re_s.LE.Re_s_thresh) THEN                 
        z0u = c_z0u_smooth*con_m/u_star     ! Smooth flow
      ELSE
        z0u = c_z0u_fetch*u_star*u_star/grav  ! Rough flow
      END IF 
      ! Roughness for scalars  
      z0q = c_z0u_fetch*MAX(Re_s, Re_s_thresh)
      z0t = c_z0t_rough_1*z0q**c_z0t_rough_3 - c_z0t_rough_2
      z0q = c_z0q_rough_1*z0q**c_z0q_rough_3 - c_z0q_rough_2
      z0t = z0u*EXP(-c_Karman/Pr_neutral*z0t)
      z0q = z0u*EXP(-c_Karman/Sc_neutral*z0q) 
      
    ELSE                               ! Ice surface 
      
      ! The Charnock parameter is not used over ice, formally set "c_z0u_fetch" to its minimum value
      c_z0u_fetch = c_z0u_rough
      
      ! Threshold value of friction velocity
      u_star_thresh = c_z0u_smooth*con_m/z0u_ice_rough
      
      ! Aerodynamic roughness
      z0u = MAX(z0u_ice_rough, c_z0u_smooth*con_m/u_star)
      
      ! Roughness Reynolds number 
      Re_s = MAX(u_star*z0u/con_m, c_accur_sf)
      
      ! Roughness for scalars  
      IF(Re_s.LE.Re_z0s_ice_t) THEN 
        z0t = c_z0t_ice_b0t + c_z0t_ice_b1t*LOG(Re_s)
        z0t = MIN(z0t, c_z0t_ice_b0s)
        z0q = c_z0q_ice_b0t + c_z0q_ice_b1t*LOG(Re_s)
        z0q = MIN(z0q, c_z0q_ice_b0s)
      ELSE 
        z0t = c_z0t_ice_b0r + c_z0t_ice_b1r*LOG(Re_s) + c_z0t_ice_b2r*LOG(Re_s)**2_wp
        z0q = c_z0q_ice_b0r + c_z0q_ice_b1r*LOG(Re_s) + c_z0q_ice_b2r*LOG(Re_s)**2_wp
      END IF
      z0t = z0u*EXP(z0t)
      z0q = z0u*EXP(z0q)
    
    END IF 

  END SUBROUTINE calc_roughness
            

      
#endif
  END MODULE mo_flake_SfcFlx_process