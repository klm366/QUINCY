!>
!! @brief Contains the routines for the lentic processes based on FLake 
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
!!  The main driving routine of the lake model FLake 
!!  where computations are performed.
!!  Advances the surface temperature
!!  and other FLake variables one time step.
!!  At the moment, the Euler explicit scheme is used.

MODULE mo_lentic_process
#ifndef __NO_JSBACH__
  
    USE mo_kind,      ONLY: wp
  
    IMPLICIT NONE
    PRIVATE
    PUBLIC ::                           &
    !  &   calc_radflux,                  &
    !  &   calc_ice_water_heatfluxes,     &
    !  &   calc_buoypar,                  &
    !  &   calc_lake_ice,                 &
    !  &   calc_lake_temp_stratification, &
     &   calc_surface_temperature

    CHARACTER(len=*), PARAMETER :: modname = 'mo_lentic_process'
    
  CONTAINS
  
    ! PURE ELEMENTAL SUBROUTINE calc_radflux ( &
    !   & depth_w,        & ! in
    !   & albedo_water,   & ! in
    !   & albedo_ice,     & ! in
    !   & albedo_snow,    & ! in
    !   & h_ice,          & ! in 
    !   & h_ML,           & ! in 
    !   & h_snow,         & ! in
    !   & I_atm,          & ! in
    !   & opticpar_water, & ! in 
    !   & opticpar_ice,   & ! in 
    !   & opticpar_snow,  & ! in 
    !   & I_bot,          & ! out
    !   & I_h,            & ! out
    !   & I_ice,          & ! out
    !   & I_intm_0_h,     & ! out
    !   & I_intm_h_D,     & ! out
    !   & I_snow,         & ! out (but maybe not needed) 
    !   & I_w             & ! out
    !   )       

    !   ! Description:
    !   !
    !   !  Computes the radiation fluxes 
    !   !  at the snow-ice, ice-water, air-water, 
    !   !  mixed layer-thermocline and water column-bottom sediment interfaces,
    !   !  the mean radiation flux over the mixed layer,
    !   !  and the mean radiation flux over the thermocline.

    !   USE mo_lentic_constants , ONLY : & 
    !     h_Snow_min            , & ! Minimum snow thickness [m]
    !     h_Ice_min             , & ! Minimum ice thickness [m]
    !     h_ML_min                  ! Minimum mixed-layer depth [m]

    !   !$claw model-data
    !   REAL (KIND = wp), INTENT(IN) ::   &
    !     depth_w       , & ! The lake depth [m]
    !     albedo_water  , & ! Albedo of the water surface 
    !     albedo_ice    , & ! Albedo of the ice surface
    !     albedo_snow   , & ! Albedo of the snow surface
    !     h_ice         , & ! Ice thickness [m]
    !     h_ML          , & ! Thickness of the mixed-layer [m] 
    !     h_snow_p      , & ! Snow thickness [m]
    !     I_atm           & ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}],
    !                       ! i.e. the incident radiation flux with no regard for the surface albedo

    !   TYPE (opticpar_medium), INTENT(IN) :: & 
    !     opticpar_water , & ! Optical characteristics of water
    !     opticpar_ice   , & ! Optical characteristics of ice
    !     opticpar_snow    & ! Optical characteristics of snow 
      
    !   REAL(wp),  INTENT(out) :: &
    !     I_bot,      & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
    !     I_h,        & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
    !     I_ice,      & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
    !     I_intm_0_h, & ! Mean radiation flux over the mixed layer [W m^{-1}]
    !     I_intm_h_D, & ! Mean radiation flux over the thermocline [W m^{-1}]
    !     I_snow,     & ! Radiation flux through the air-snow interface [W m^{-2}] (but maybe not needed) 
    !     I_w,        & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
    !   !$claw end model-data


    !   !  Local variables of type INTEGER
    !   INTEGER (KIND = iintegers) :: & ! Help variable(s)
    !   i                               ! DO loop index

    !   ! ---------------------------
    !   ! Go

    !   IF(h_ice.GE.h_Ice_min) THEN            ! Ice exists
    !     IF(h_snow_p.GE.h_Snow_min) THEN        ! There is snow above the ice
    !       I_snow = I_atm*(1._wp-albedo_snow) 
    !       I_bot = 0._wp
    !       DO i=1, opticpar_snow%nband_optic
    !         I_bot = I_bot +                    & 
    !         opticpar_snow%frac_optic(i)*EXP(-opticpar_snow%extincoef_optic(i)*h_snow_p) 
    !       END DO 
    !       I_ice  = I_snow*I_bot
    !     ELSE                                           ! No snow above the ice 
    !       I_snow = I_atm  
    !       I_ice  = I_atm*(1._wp-albedo_ice)
    !     END IF 
    !     I_bot = 0._wp
    !     DO i=1, opticpar_ice%nband_optic
    !       I_bot = I_bot +                      & 
    !       opticpar_ice%frac_optic(i)*EXP(-opticpar_ice%extincoef_optic(i)*h_ice) 
    !     END DO 
    !     I_w      = I_ice*I_bot
    !   ELSE                                             ! No ice-snow cover
    !     I_snow   = I_atm  
    !     I_ice    = I_atm
    !     I_w      = I_atm*(1._wp-albedo_water)
    !   END IF 

    !   IF(h_ML.GE.h_ML_min) THEN           ! Radiation flux at the bottom of the mixed layer
    !     I_bot = 0._wp
    !     DO i=1, opticpar_water%nband_optic
    !       I_bot = I_bot +            & 
    !       opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*h_ML) 
    !     END DO 
    !     I_h = I_w*I_bot
    !   ELSE                                          ! Mixed-layer depth is less then a minimum value
    !     I_h = I_w
    !   END IF

    !   I_bot = 0._wp                         ! Radiation flux at the lake bottom
    !   DO i=1, opticpar_water%nband_optic
    !     I_bot = I_bot +              & 
    !     opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*depth_w) 
    !   END DO 
    !   I_bot = I_w*I_bot

    !   IF(h_ML.GE.h_ML_min) THEN           ! Integral-mean radiation flux over the mixed layer
    !     I_intm_0_h = 0._wp
    !     DO i=1, opticpar_water%nband_optic
    !       I_intm_0_h = I_intm_0_h +                                &
    !       opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
    !       (1._wp - EXP(-opticpar_water%extincoef_optic(i)*h_ML))
    !     END DO 
    !     I_intm_0_h = I_w*I_intm_0_h/h_ML
    !   ELSE
    !     I_intm_0_h = I_h
    !   END IF

    !   IF(h_ML.LE.depth_w-h_ML_min) THEN   ! Integral-mean radiation flux over the thermocline
    !     I_intm_h_D = 0._wp 
    !     DO i=1, opticpar_water%nband_optic
    !       I_intm_h_D = I_intm_h_D +                                &
    !       opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
    !       ( EXP(-opticpar_water%extincoef_optic(i)*h_ML)             &
    !       - EXP(-opticpar_water%extincoef_optic(i)*depth_w) )
    !     END DO 
    !     I_intm_h_D = I_w*I_intm_h_D/(depth_w-h_ML)
    !   ELSE
    !     I_intm_h_D = I_h
    !   END IF

    ! END SUBROUTINE flake_radflux

    ! PURE ELEMENTAL SUBROUTINE calc_buoypar (&
    !   & T_water,       & ! in
    !   & buoypar        & ! out
    !   )

    !   ! Description:
    !   !  Computes the buoyancy parameter,
    !   !  using a quadratic equation of state for the fresh-water.
      
    !   USE mo_lentic_constants,  ONLY: &
    !     & tpl_grav                  , & ! Acceleration due to gravity [m s^{-2}]
    !     & tpl_T_r                   , & ! Temperature of maximum density of fresh water [K]
    !     & tpl_a_T                       ! Constant in the fresh-water equation of state [K^{-2}]
      

    !   !$claw model-data
    !   REAL(wp), INTENT(in)  :: T_water   ! Water temperature [K]
    !   REAL(wp), INTENT(out) :: buoypar   ! Buoyancy parameter [m s^{-2} K^{-1}]
    !   !$claw end model-data

    !   ! ---------------------------
    !   ! Go

    !   buoypar = tpl_grav*tpl_a_T*(T_water-tpl_T_r)
      
    ! END SUBROUTINE calc_buoypar

    ! PURE ELEMENTAL SUBROUTINE calc_ice_water_heatfluxes( &
    !   & h_ice,      & ! in
    !   & h_ML,       & ! in
    !   & T_bot,      & ! in
    !   & T_wML,      & ! in
    !   & C_T,        & ! in
    !   & depth_w,    & ! in
    !   & I_w,        & ! in
    !   & I_h,        & ! in
    !   & I_intm_0_h, & ! in
    !   & depth_w,    & ! in
    !   & I_bot,      & ! inout
    !   & Q_w,        & ! out
    !   & Q_star,     & ! out
    !   & Phi_T_pr0   & ! out
    !   )

    !   USE mo_lentic_constants,  ONLY: &
    !     & h_Ice_min,     &
    !     & h_ML_min,      &
    !     & tpl_kappa_w,   &
    !     & Phi_T_pr0_1,   &
    !     & Phi_T_pr0_2,   &
    !     & albedo_bot    
        
        
    !   !$claw model-data
    !   REAL(wp), INTENT(in) :: h_ice,      & ! Ice thickness [m]
    !                         & h_ML,       & ! Thickness of the mixed-layer [m] 
    !                         & T_bot,      & ! Temperature at the water-bottom sediment interface [K] 
    !                         & T_wML,      & ! Mixed-layer temperature [K] 
    !                         & C_T,        & ! Shape factor (thermocline)
    !                         & I_w,        & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
    !                         & I_h,        & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
    !                         & I_intm_0_h, & ! Mean radiation flux over the mixed layer [W m^{-1}]
    !                         & depth_w       ! Waterbody depth [m]

    !   REAL(wp),  INTENT(inout) :: I_bot  ! Radiation flux through the water-bottom sediment interface [W m^{-2}]

    !   REAL(wp),  INTENT(out) :: Q_w,       & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
    !                           & Q_star,    & ! A generalized heat flux scale [W m^{-2}]
    !                           & Phi_T_pr0    ! d\Phi_T(0)/d\zeta   (thermocline)
    !   !$claw end model-data

    !   ! ---------------------------
    !   ! Go

    !   ! Compute heat flux through the ice-water interface
    !   IF(h_ice.GE.h_Ice_min) THEN                   ! Ice exists 
    !     IF(h_ML.LE.h_ML_min) THEN                   ! Mixed-layer depth is zero, compute flux 
    !       Q_w = -tpl_kappa_w*(T_bot-T_wML)/depth_w  ! Flux with linear T(z) 
    !       Phi_T_pr0 = Phi_T_pr0_1*C_T-Phi_T_pr0_2   ! d\Phi(0)/d\zeta (thermocline)
    !       Q_w = Q_w*MAX(Phi_T_pr0, 1._wp)           ! Account for an increased d\Phi(0)/d\zeta 
    !     ELSE                    
    !       Q_w = 0._wp                               ! Mixed-layer depth is greater than zero, set flux to zero
    !     END IF   
    !   END IF   
        
    !   ! A generalized heat flux scale 
    !   I_bot = (1._wp - albedo_bot)*I_bot
    !   Q_star = Q_w + I_w + I_h - 2._wp*I_intm_0_h         
    ! END SUBROUTINE calc_ice_water_heatfluxes

    ! PURE ELEMENTAL SUBROUTINE calc_lake_ice( &
    !   & I_ice,    & ! in
    !   & I_w,      & ! in 
    !   & T_wML,    & ! in 
    !   & Q_ice,    & ! in
    !   & Q_w,      & ! in
    !   & del_time, & ! in
    !   & h_ice,    & ! inout
    !   & T_ice     & ! inout
    !   & l_ice_create & ! out
    !   )
    
    !   USE mo_lentic_constants,  ONLY: &
    !     & C_I_lin,       &
    !     & C_I_MR,        &
    !     & h_Ice_min,     &
    !     & tpl_c_I,       &
    !     & tpl_kappa_I,   &
    !     & tpl_L_f,       &
    !     & tpl_T_f,       &
    !     & tpl_rho_I,     &
    !     & c_small,       &
    !     & Phi_I_pr0_lin, &
    !     & Phi_I_pr1_lin, &
    !     & Phi_I_ast_MR,  &
    !     & H_Ice_max     
        
    !   !$claw model-data      
    !   REAL(wp), INTENT(in) ::del_time, & ! The model time step [s]
    !                         &   h_ice, & ! Ice thickness [m]
    !                         &   T_ice, & ! Temperature at the snow-ice or air-ice interface [K] 
    !                         &   I_ice, & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
    !                         &   I_w,   & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
    !                         &   T_wML, & ! Mixed-layer temperature [K] 
    !                         &   Q_ice, & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]   
    !                         &   Q_w      ! Heat flux through the ice-water or air-water interface [W m^{-2}]
    
    !   LOGICAL, INTENT(out) :: l_ice_create ! Switch, .TRUE. = ice does not exist but should be created
    
    !   ! Local variables
    !   LOGICAL :: l_ice_meltabove ! Switch, .TRUE. = snow/ice melting from above takes place
    
    !   REAL(wp) :: d_h_ice_dt,     & ! Time derivative of h_ice [m s^{-1}]
    !             & C_I,            & ! Shape factor (ice)
    !             & h_ice_threshold,& ! If h_ice<h_ice_threshold, use quasi-equilibrium ice model 
    !             & h_ice_p,        & ! Ice thickness in previous timestep [m]
    !             & Phi_I_pr0,      & ! d\Phi_I(0)/d\zeta_I (ice)
    !             & Phi_I_pr1       & ! d\Phi_I(1)/d\zeta_I (ice)
    !             & flk_str_1,      & ! Help storage variable
    !             & flk_str_2,      & ! Help storage variable
    !             & R_TI_icesnow,   & ! Dimensionless ratio, used to store intermediate results
    !             & R_Tstar_icesnow,& ! Dimensionless ratio, used to store intermediate results
    
    
    !   ! Initialization
    !   l_ice_create    = .FALSE.  
    !   l_ice_meltabove = .FALSE.  
    
    !   d_h_ice_dt   = 0._wp 
    !   h_ice_p      = h_ice 

    !   ! ---------------------------
    !   ! Go
    
    !   IF(h_ice_p.LT.h_Ice_min) THEN   
    !   ! Ice does not exist 
      
    !     l_ice_create = T_wML.LE.(tpl_T_f+c_small).AND.Q_w.LT.0._wp
    
    !     IF(l_ice_create) THEN                            
    !     ! Ice does not exist but should be created
    !       d_h_ice_dt = -Q_w/tpl_rho_I/tpl_L_f                                  
    !       h_ice = h_ice_p + d_h_ice_dt*del_time                    ! Advance h_ice 
    !       T_ice = tpl_T_f + h_ice*Q_w/tpl_kappa_I/Phi_I_pr0_lin    ! Ice temperature
    !       Phi_I_pr1 = Phi_I_pr1_lin                                    & 
    !                     + Phi_I_ast_MR*MIN(1._wp, h_ice/H_Ice_max) ! d\Phi_I(1)/d\zeta_I (ice)
    !     END IF
    !   ELSE                                
    !   ! Ice exists
    !     IF(T_ice.GE.(tpl_T_f-c_small)) THEN  
    !     ! T_sfc = T_f, check for melting from above     
    !       flk_str_1 = Q_ice + I_ice - I_w - Q_w  ! Atmospheric forcing + heating from the water
          
    !       IF(flk_str_1.GE.0._wp) THEN  
    !       ! Melting of ice from above
    !         l_ice_meltabove = .TRUE.
    !         d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
    !         h_ice  = h_ice_p  + d_h_ice_dt *del_time  ! Advance h_ice
    !         T_ice  = tpl_T_f                          ! Set T_ice to the freezing point
    !       END IF
        
    !     END IF 
        
    !     IF(.NOT.l_ice_meltabove) THEN                 
    !     ! No melting from above
          
    !       flk_str_1 = 0._wp
        
    !       Phi_I_pr0 = h_ice_p/H_Ice_max                          ! h_ice relative to its maximum value
    !       C_I = C_I_lin - C_I_MR*(1._wp+Phi_I_ast_MR)*Phi_I_pr0  ! Shape factor (ice)
    !       Phi_I_pr1 = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr0     ! d\Phi_I(1)/d\zeta_I (ice)
    !       Phi_I_pr0 = Phi_I_pr0_lin - Phi_I_pr0                  ! d\Phi_I(0)/d\zeta_I (ice)
          
    !       h_ice_threshold = MAX(1._wp, 2._wp*C_I*tpl_c_I*(tpl_T_f-T_ice)/tpl_L_f)
    !       h_ice_threshold = Phi_I_pr0/C_I*tpl_kappa_I/tpl_rho_I/tpl_c_I*h_ice_threshold
    !       h_ice_threshold = SQRT(h_ice_threshold*del_time)       ! Threshold value of h_ice
    !       h_ice_threshold = MIN(0.._wp*H_Ice_max, MAX(h_ice_threshold, h_Ice_min))
    !                                                             ! h_ice(threshold) < 0.9*H_Ice_max
          
    !       IF(h_ice_p.LT.h_ice_threshold) THEN  
    !       ! Use a quasi-equilibrium ice model
    !         flk_str_1 = Q_ice + I_ice - I_w
    !         d_h_ice_dt = -(flk_str_1-Q_w)/tpl_L_f/tpl_rho_I
    !         h_ice = h_ice_p + d_h_ice_dt *del_time                   ! Advance h_ice
    !         T_ice = tpl_T_f + h_ice*flk_str_1/tpl_kappa_I/Phi_I_pr0  ! Ice temperature    
    !       ELSE                                    
    !       ! Use a complete ice model     
    !         d_h_ice_dt = tpl_kappa_I*(tpl_T_f-T_ice)/h_ice_p*Phi_I_pr0
    !         d_h_ice_dt = (Q_w+d_h_ice_dt)/tpl_L_f/tpl_rho_I
    !         h_ice = h_ice_p  + d_h_ice_dt*del_time          ! Advance h_ice
            
    !         R_TI_icesnow = tpl_c_I*(tpl_T_f-T_ice)/tpl_L_f  ! Dimensionless parameter
    !         R_Tstar_icesnow = (1._wp - C_I)*R_TI_icesnow    ! Dimensionless parameter
    !         flk_str_2  = Q_ice+I_ice-I_w                    ! Atmospheric fluxes
    !         flk_str_1  = C_I*h_ice_p
    !         d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-T_ice)/h_ice_p*Phi_I_pr0       &
    !                   * (1._wp-R_Tstar_icesnow)             ! Add flux due to heat conduction
    !         d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*Q_w   ! Add flux from water to ice
    !         d_T_ice_dt = d_T_ice_dt + flk_str_2             ! Add atmospheric fluxes
    !         d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I       ! Total forcing
    !         d_T_ice_dt = d_T_ice_dt/flk_str_1               ! dT_ice/dt 
    !         T_ice = T_ice + d_T_ice_dt*del_time             ! Advance T_ice
    !       END IF
          
    !       Phi_I_pr1 = MIN(1._wp, h_ice/H_Ice_max)             ! h_ice relative to its maximum value
    !       Phi_I_pr1 = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr1  ! d\Phi_I(1)/d\zeta_I (ice)
          
    !     END IF 
      
    !   END IF    
      
    !   ! Security, limit h_ice by its maximum value
    !   h_ice = MIN(h_ice, H_Ice_max)      
      
    !   ! Security, limit the ice and snow temperatures by the freezing point 
    !   T_ice =  MIN(T_ice,  tpl_T_f)    
      
    !   !_tmp
    !   ! Security, avoid too low values (these constraints are used for debugging purposes)
    !   T_ice =  MAX(T_ice,  73.1._wp)    
      
    !   ! Remove too thin ice and/or snow
    !   IF(h_ice.LT.h_Ice_min)  THEN  ! Check ice
    !     h_ice = 0._wp               ! Ice is too thin, remove it, and
    !     T_ice = tpl_T_f             ! set T_ice to the freezing point.
      
    !     l_ice_create = .FALSE.      ! "Exotic" case, ice has been created but proved to be too thin
    !   END IF
    ! END SUBROUTINE calc_lake_ice

    ! PURE ELEMENTAL SUBROUTINE calc_lake_temp_stratification( &
    !   & del_time,             & ! in
    !   & depth_w,              & ! in
    !   & extincoef_water_typ,  & ! in
    !   & l_ice_create,         & ! in
    !   & h_ice,                & ! in
    !   & I_bot,                & ! in
    !   & I_intm_0_h,           & ! in
    !   & I_w,                  & ! in
    !   & par_Coriolis,         & ! in
    !   & Q_w,                  & ! in
    !   & Q_bot,                & ! in
    !   & Q_star,               & ! in
    !   & u_star_w,             & ! in 
    !   & C_T,                  & ! inout
    !   & T_bot,                & ! inout
    !   & T_mnw,                & ! inout
    !   & T_wML,                & ! inout
    !   & h_ML,                 & ! inout
    !   & w_star_sfc              ! out
    !   )

    !   USE mo_lentic_constants,  ONLY: &
    !     & c_cbl_1,      &
    !     & c_cbl_2,      &
    !     & c_sbl_ZM_i,   &
    !     & c_sbl_ZM_n,   &
    !     & c_sbl_ZM_s,   &
    !     & c_small,      &
    !     & c_relax_C,    &
    !     & c_relax_h,    & 
    !     & C_T_min,      &
    !     & C_T_max,      &
    !     & C_TT_1,       &
    !     & C_TT_2,       &
    !     & h_ML_min,     &
    !     & h_ML_max,     & 
    !     & h_Ice_min,    &
    !     & tpl_c_w,      &
    !     & tpl_rho_w_r,  &
    !     & tpl_T_f,      &
    !     & tpl_T_r,      &
    !     & u_star_min

    !   !$claw model-data
    !   REAL(wp), INTENT(in) :: del_time,             & ! The model time step [s]
    !                         & depth_w,              & ! The waterbody depth [m]
    !                         & extincoef_water_typ,  & ! "Typical" extinction coefficient of the lake water [m^{-1}], used to compute the equilibrium CBL depth
    !                         & l_ice_create,         & ! Switch, .TRUE. = ice does not exist but should be created
    !                         & h_ice,                & ! Ice thickness [m]
    !                         & I_bot,                & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
    !                         & I_intm_0_h,           & ! Mean radiation flux over the mixed layer [W m^{-1}]
    !                         & I_w,                  & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
    !                         & par_Coriolis,         & ! The Coriolis parameter [s^{-1}]
    !                         & Q_w,                  & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
    !                         & Q_bot,                & ! Heat flux through the water-bottom sediment interface [W m^{-2}]
    !                         & Q_star,               & ! A generalized heat flux scale [W m^{-2}]
    !                         & u_star_w                ! Friction velocity in the surface layer of lake water [m s^{-1}]

    !   REAL(wp),  INTENT(inout) :: C_T,                  & ! Shape factor (thermocline)
    !                             & T_bot,                & ! Temperature at the water-bottom sediment interface [K] 
    !                             & T_mnw,                & ! Mean temperature of the water column [K] 
    !                             & T_wML,                & ! Mixed-layer temperature [K] 
    !                             & h_ML                    ! Thickness of the mixed-layer [m] 

    !   REAL(wp),  INTENT(out) :: w_star_sfc ! Convective velocity scale, using a generalized heat flux scale [m s^{-1}]
    !   !$claw end model-data

    !   ! Local variables
    !   REAL(wp) ::
    !     & C_T_p,              & ! Shape factor in previous timestep (thermocline)
    !     & T_bot_p,            & ! Temperature at the water-bottom sediment interface in previous timestep [K] 
    !     & T_wML_p,            & ! Mixed-layer temperature in previous timestep [K] 
    !     & h_ML_p,             & ! Thickness of the mixed-layer in previous timestep [m] 
    !     & C_Q,                & ! Shape factor with respect to the heat flux (thermocline)
    !     & C_TT,               & ! Dimensionless parameter (thermocline)
    !     & d_T_mnw_dt,         & ! Time derivative of T_mnw [K s^{-1}] 
    !     & d_T_bot_dt,         & ! Time derivative of T_bot [K s^{-1}] 
    !     & d_h_ML_dt,          & ! Time derivative of h_ML [m s^{-1}]
    !     & flk_str_1,          & ! Help storage variable
    !     & flk_str_2,          & ! Help storage variable
    !     & buoypar,            & ! Buoyancy parameter [m s^{-2} K^{-1}]
    !     & conv_equil_h_scale, & ! The equilibrium CBL depth scale [m]
    !     & N_T_mean,           & ! The mean buoyancy frequency in the thermocline [s^{-1}] 
    !     & d_C_T_dt,           & ! Time derivative of C_T [s^{-1}]
    !     & w_star_sfc,         & ! Convective velocity scale, using a generalized heat flux scale [m s^{-1}]
    !     & R_H_icesnow,        & ! Dimensionless ratio, used to store intermediate results
    !     & R_rho_c_icesnow,    & ! Dimensionless ratio, used to store intermediate results
    !     & R_TI_icesnow,       & ! Dimensionless ratio, used to store intermediate results
    !     & R_Tstar_icesnow,    & ! Dimensionless ratio, used to store intermediate results
    !     & ZM_h_scale

    !   ! Initialization
    !   C_T_p   = C_T
    !   T_bot_p = T_bot
    !   T_wML_p = T_wML
    !   h_ML_p  = h_ML
    !   d_T_mnw_dt   = 0._wp 
    !   d_T_bot_dt   = 0._wp 
    !   d_h_ML_dt    = 0._wp 
    !   d_C_T_dt     = 0._wp 

    !   ! ---------------------------
    !   ! Go

    !   !  Compute the mean temperature of the water column.
    !   IF(l_ice_create) Q_w = 0._wp     ! Ice has just been created, set Q_w to zero
    !   d_T_mnw_dt = (Q_w - Q_bot + I_w - I_bot)/tpl_rho_w_r/tpl_c_w/depth_w

    !   T_mnw = T_mnw + d_T_mnw_dt*del_time   ! Advance T_mnw
    !   T_mnw = MAX(T_mnw, tpl_T_f)           ! Limit T_mnw by the freezing point 
    
    
    !   !------------------------------------------------------------------------------
    !   !  Compute the mixed-layer depth, the mixed-layer temperature, 
    !   !  the bottom temperature and the shape factor
    !   !  with respect to the temperature profile in the thermocline. 
    !   !  Different formulations are used, depending on the regime of mixing. 
    !   !------------------------------------------------------------------------------
    
    !   IF(h_ice.GE.h_Ice_min) THEN    
    !     ! Ice exists
    
    !     T_mnw = MIN(T_mnw, tpl_T_r) ! Limit the mean temperature under the ice by T_r 
    !     T_wML = tpl_T_f                   ! The mixed-layer temperature is equal to the freezing point 
    
    !     IF(l_ice_create) THEN                  ! Ice has just been created 
    !       IF(h_ML_p.GE.depth_w-h_ML_min) THEN    ! h_ML=D when ice is created 
    !         h_ML = 0._wp                 ! Set h_ML to zero 
    !         C_T = C_T_min                    ! Set C_T to its minimum value 
    !       ELSE                                          ! h_ML<D when ice is created 
    !         h_ML = h_ML_p                ! h_ML remains unchanged 
    !         ! C_T = C_T_p, C_T (thermocline) remains unchanged 
    !       END IF 
    !         T_bot = T_wML - (T_wML-T_mnw)/C_T/(1._wp-h_ML/depth_w)
    !                             ! Update the bottom temperature 
        
    !     ELSE IF(T_bot_p.LT.tpl_T_r) THEN   ! Ice exists and T_bot < T_r, molecular heat transfer 
    !       h_ML = h_ML_p                  ! h_ML remains unchanged 
    !       ! C_T = C_T_p, C_T (thermocline) remains unchanged 
    !       T_bot = T_wML - (T_wML-T_mnw)/C_T/(1._wp-h_ML/depth_w)
    !                                 ! Update the bottom temperature 
            
    !     ELSE                                   ! Ice exists and T_bot = T_r, convection due to bottom heating 
    !       T_bot = tpl_T_r                      ! T_bot is equal to the temperature of maximum density 
    !       IF(h_ML_p.GE.c_small) THEN   ! h_ML > 0 
    !         ! C_T = C_T_p, C_T (thermocline) remains unchanged 
    !         h_ML = depth_w*(1._wp-(T_wML-T_mnw)/(T_wML-T_bot)/C_T)
    !         h_ML = MAX(h_ML, 0._wp)   ! Update the mixed-layer depth  
    !       ELSE                                 ! h_ML = 0 
    !         h_ML = h_ML_p                   ! h_ML remains unchanged 
    !         C_T = (T_wML-T_mnw)/(T_wML-T_bot) 
    !         C_T = MIN(C_T_max, MAX(C_T, C_T_min)) ! Update the shape factor (thermocline)  
    !       END IF 
    !     END IF 
      
    !     T_bot = MIN(T_bot, tpl_T_r)    ! Security, limit the bottom temperature by T_r 
    !   ELSE
    !     ! Open water
    
    !     ! Generalised buoyancy flux scale and convective velocity scale
    !     buoypar = calc_buoypar(T_wML_p)
    !     flk_str_1 = buoypar*Q_star/tpl_rho_w_r/tpl_c_w                    
    !     IF(flk_str_1.LT.0._wp) THEN       
    !       w_star_sfc = (-flk_str_1*h_ML_p)**(1._wp/3._wp)  ! Convection     
    !     ELSE 
    !       w_star_sfc = 0._wp                                       ! Neutral or stable stratification
    !     END IF 
        
    !     !_dm
    !     ! The equilibrium depth of the CBL due to surface cooling with the volumetric heating
    !     ! is not computed as a solution to the transcendental equation.
    !     ! Instead, an algebraic formula is used
    !     ! that interpolates between the two asymptotic limits.
    !     !_dm
    !     conv_equil_h_scale = -Q_w/MAX(I_w, c_small)
    !     IF(conv_equil_h_scale.GT.0._wp .AND. conv_equil_h_scale.LT.1._wp  &
    !     .AND. T_wML_p.GT.tpl_T_r) THEN   ! The equilibrium CBL depth scale is only used above T_r
    !       conv_equil_h_scale = SQRT(6._wp*conv_equil_h_scale)                 &
    !         + 2._wp*conv_equil_h_scale/(1._wp-conv_equil_h_scale)
    !       conv_equil_h_scale = MIN(depth_w, conv_equil_h_scale/extincoef_water_typ)
    !     ELSE
    !       conv_equil_h_scale = 0._wp       ! Set the equilibrium CBL depth to zero
    !     END IF
    
    !     ! Mean buoyancy frequency in the thermocline
    !     buoypar = calc_buoypar(0.5_wp*(T_wML_p+T_bot_p))
    !     N_T_mean = buoypar*(T_wML_p-T_bot_p)
    !     IF(h_ML_p.LE.depth_w-h_ML_min) THEN
    !       N_T_mean = SQRT(N_T_mean/(depth_w-h_ML_p))  ! Compute N                   
    !     ELSE 
    !       N_T_mean = 0._wp                            ! h_ML=D, set N to zero
    !     END IF 
    
    !     ! The rate of change of C_T
    !     d_C_T_dt = MAX(w_star_sfc, u_star_w, u_star_min)**2_iintegers
    !     d_C_T_dt = N_T_mean*(depth_w-h_ML_p)**2_iintegers       &
    !     / c_relax_C/d_C_T_dt                               ! Relaxation time scale for C_T
    !     d_C_T_dt = (C_T_max-C_T_min)/MAX(d_C_T_dt, c_small)     ! Rate-of-change of C_T 
        
    !     ! Compute the shape factor and the mixed-layer depth, 
    !     ! using different formulations for convection and wind mixing
        
    !     C_TT = C_TT_1*C_T_p-C_TT_2         ! C_TT, using C_T at the previous time step
    !     C_Q = 2._wp*C_TT/C_T_p     ! C_Q using C_T at the previous time step
        
    !     IF(flk_str_1.LT.0._wp) THEN  ! Convective mixing 
    
    !       C_T = C_T_p + d_C_T_dt*del_time                        ! Update C_T, assuming dh_ML/dt>0
    !       C_T = MIN(C_T_max, MAX(C_T, C_T_min))                ! Limit C_T 
    !       d_C_T_dt = (C_T-C_T_p)/del_time                        ! Re-compute dC_T/dt
          
    !       IF(h_ML_p.LE.depth_w-h_ML_min) THEN       ! Compute dh_ML/dt
    !         IF(h_ML_p.LE.h_ML_min) THEN    ! Use a reduced entrainment equation (spin-up)
    !           d_h_ML_dt = c_cbl_1/c_cbl_2*MAX(w_star_sfc, c_small)
              
    !         ELSE                                   ! Use a complete entrainment equation 
    !           R_H_icesnow     = depth_w/h_ML_p
    !           R_rho_c_icesnow = R_H_icesnow-1._wp
    !           R_TI_icesnow    = C_T_p/C_TT
    !           R_Tstar_icesnow = (R_TI_icesnow/2._wp-1._wp)*R_rho_c_icesnow + 1._wp
    !           d_h_ML_dt = -Q_star*(R_Tstar_icesnow*(1._wp+c_cbl_1)-1._wp) - Q_bot
    !           d_h_ML_dt = d_h_ML_dt/tpl_rho_w_r/tpl_c_w                        ! Q_* and Q_b flux terms
    !           flk_str_2 = (depth_w-h_ML_p)*(T_wML_p-T_bot_p)*C_TT_2/C_TT*d_C_T_dt 
    !           d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add dC_T/dt term
    !           flk_str_2 = I_bot + (R_TI_icesnow-1._wp)*I_h - R_TI_icesnow*I_intm_h_D
    !           flk_str_2 = flk_str_2 + (R_TI_icesnow-2._wp)*R_rho_c_icesnow*(I_h-I_intm_0_h)
    !           flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w
    !           d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add radiation terms
    !           flk_str_2 = -c_cbl_2*R_Tstar_icesnow*Q_star/tpl_rho_w_r/tpl_c_w/MAX(w_star_sfc, c_small)
    !           flk_str_2 = flk_str_2 + C_T_p*(T_wML_p-T_bot_p)
    !           d_h_ML_dt = d_h_ML_dt/flk_str_2                                   ! dh_ML/dt = r.h.s.
    !         END IF 
    !         !_dm
    !         ! Notice that dh_ML/dt may appear to be negative  
    !         ! (e.g. due to buoyancy loss to bottom sediments and/or
    !         ! the effect of volumetric radiation heating),
    !         ! although a negative generalized buoyancy flux scale indicates 
    !         ! that the equilibrium CBL depth has not yet been reached
    !         ! and convective deepening of the mixed layer should take place.
    !         ! Physically, this situation reflects an approximate character of the lake model.
    !         ! Using the self-similar temperature profile in the thermocline, 
    !         ! there is always communication between the mixed layer, the thermocline 
    !         ! and the lake bottom. As a result, the rate of change of the CBL depth
    !         ! is always dependent on the bottom heat flux and the radiation heating of the thermocline.
    !         ! In reality, convective mixed-layer deepening may be completely decoupled
    !         ! from the processes underneath. In order to account for this fact,
    !         ! the rate of CBL deepening is set to a small value
    !         ! if dh_ML/dt proves to be negative.
    !         ! This is "double insurance" however, 
    !         ! as a negative dh_ML/dt is encountered very rarely.
    !         !_dm
    
    
    !         d_h_ML_dt = MAX(d_h_ML_dt, c_small)    
    !         h_ML = h_ML_p + d_h_ML_dt*del_time                       ! Update h_ML 
    !         h_ML = MAX(h_ML_min, MIN(h_ML, depth_w))           ! Security, limit h_ML
    !       ELSE                                              ! Mixing down to the lake bottom
    !         h_ML = depth_w
    !       END IF
    
    !     ELSE  
    !       ! Wind mixing
    
    !       d_h_ML_dt = MAX(u_star_w, u_star_min)                        ! The surface friction velocity
    !       ZM_h_scale = (ABS(par_Coriolis)/c_sbl_ZM_n + N_T_mean/c_sbl_ZM_i)*d_h_ML_dt**2_iintegers
    !       ZM_h_scale = ZM_h_scale + flk_str_1/c_sbl_ZM_s
    !       ZM_h_scale = MAX(ZM_h_scale, c_small)
    !       ZM_h_scale = d_h_ML_dt**3_iintegers/ZM_h_scale 
    !       ZM_h_scale = MAX(h_ML_min, MIN(ZM_h_scale, h_ML_max))        ! The ZM96 SBL depth scale 
    !       ZM_h_scale = MAX(ZM_h_scale, conv_equil_h_scale)                     ! Equilibrium mixed-layer depth 
          
    !       !_dm 
    !       ! In order to avoid numerical discretization problems,
    !       ! an analytical solution to the evolution equation 
    !       ! for the wind-mixed layer depth is used.
    !       ! That is, an exponential relaxation formula is applied
    !       ! over the time interval equal to the model time step.
    !       !_dm 
          
    !       d_h_ML_dt = c_relax_h*d_h_ML_dt/ZM_h_scale*del_time
    !       h_ML = ZM_h_scale - (ZM_h_scale-h_ML_p)*EXP(-d_h_ML_dt)    ! Update h_ML 
    !       h_ML = MAX(h_ML_min, MIN(h_ML, depth_w))             ! Limit h_ML 
    !       d_h_ML_dt = (h_ML-h_ML_p)/del_time                         ! Re-compute dh_ML/dt
          
    !       IF(h_ML.LE.h_ML_p)           &
    !         d_C_T_dt = -d_C_T_dt                 ! Mixed-layer retreat or stationary state, dC_T/dt<0
    !       C_T = C_T_p + d_C_T_dt*del_time                            ! Update C_T
    !       C_T = MIN(C_T_max, MAX(C_T, C_T_min))                    ! Limit C_T 
    !       d_C_T_dt = (C_T-C_T_p)/del_time                            ! Re-compute dC_T/dt
                  
    !     END IF 
    
    !     ! Compute the time-rate-of-change of the the bottom temperature, 
    !     ! depending on the sign of dh_ML/dt 
    !     ! Update the bottom temperature and the mixed-layer temperature
        
    !     IF(h_ML.LE.depth_w-h_ML_min) THEN       ! Mixing did not reach the bottom 
        
    !       IF(h_ML.GT.h_ML_p) THEN   
    !         ! Mixed-layer deepening 
    !         R_H_icesnow     = h_ML_p/depth_w
    !         R_rho_c_icesnow = 1._wp-R_H_icesnow 
    !         R_TI_icesnow    = 0.5_wp*C_T_p*R_rho_c_icesnow+C_TT*(2._wp*R_H_icesnow-1._wp)
    !         R_Tstar_icesnow = (0.5_wp+C_TT-C_Q)/R_TI_icesnow
    !         R_TI_icesnow    = (1._wp-C_T_p*R_rho_c_icesnow)/R_TI_icesnow
            
    !         d_T_bot_dt = (Q_w-Q_bot+I_w-I_bot)/tpl_rho_w_r/tpl_c_w
    !         d_T_bot_dt = d_T_bot_dt - C_T_p*(T_wML_p-T_bot_p)*d_h_ML_dt
    !         d_T_bot_dt = d_T_bot_dt*R_Tstar_icesnow/depth_w                   ! Q+I fluxes and dh_ML/dt term
            
    !         flk_str_2 = I_intm_h_D - (1._wp-C_Q)*I_h - C_Q*I_bot
    !         flk_str_2 = flk_str_2*R_TI_icesnow/(depth_w-h_ML_p)/tpl_rho_w_r/tpl_c_w
    !         d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add radiation-flux term
            
    !         flk_str_2 = (1._wp-C_TT_2*R_TI_icesnow)/C_T_p
    !         flk_str_2 = flk_str_2*(T_wML_p-T_bot_p)*d_C_T_dt
    !         d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add dC_T/dt term    
    !       ELSE                                ! Mixed-layer retreat or stationary state
    !         d_T_bot_dt = 0._wp                                            ! dT_bot/dt=0
    !       END IF
          
    !       T_bot = T_bot_p + d_T_bot_dt*del_time                      ! Update T_bot  
    !       T_bot = MAX(T_bot, tpl_T_f)           ! Security, limit T_bot by the freezing point
    !       buoypar = calc_buoypar(T_mnw)
    !       flk_str_2 = (T_bot-tpl_T_r)*buoypar

    !       IF(flk_str_2.LT.0._wp) T_bot = tpl_T_r  ! Security, avoid T_r crossover 
    !       T_wML = C_T*(1._wp-h_ML/depth_w)
    !       T_wML = (T_mnw-T_bot*T_wML)/(1._wp-T_wML)
    !       T_wML = MAX(T_wML, tpl_T_f)           ! Security, limit T_wML by the freezing point 
        
    !     ELSE                                              
    !       ! Mixing down to the lake bottom 
    !       h_ML = depth_w
    !       T_wML = T_mnw
    !       T_bot = T_mnw
    !       C_T = C_T_min      
    !     END IF
        
    !   END IF 

    !   ! In case of unstable stratification, force mixing down to the bottom
    !   buoypar = calc_buoypar(T_mnw)
    !   flk_str_2 = (T_wML-T_bot)*buoypar
    !   IF(flk_str_2.LT.0._wp) THEN 
    !     h_ML = depth_w
    !     T_wML = T_mnw
    !     T_bot = T_mnw
    !     C_T = C_T_min
    !   END IF

    ! END SUBROUTINE calc_lake_temp_stratification


    PURE ELEMENTAL SUBROUTINE calc_surface_temperature(T_sfc )!&
      ! & h_snow, & ! in
      ! & T_snow, & ! in
      ! & h_ice,  & ! in
      ! & T_ice,  & ! in
      ! & T_wML,  & ! in
      ! & T_sfc   ) ! out

      ! USE mo_lentic_constants,  ONLY: &
      !   & h_Snow_min, &
      !   & h_Ice_min   

      !$claw model-data
      ! REAL(wp), INTENT(in)  :: h_snow  ! Snow thickness [m]
      ! REAL(wp), INTENT(in)  :: T_snow  ! Snow temperature [K]
      ! REAL(wp), INTENT(in)  :: h_ice   ! Ice thickness [m]
      ! REAL(wp), INTENT(in)  :: T_ice   ! Ice temperature [K]
      ! REAL(wp), INTENT(in)  :: T_wML   ! Water temperature [K]
      REAL(wp), INTENT(out) :: T_sfc   ! Surface temperature [K]
      !$claw end model-data

      ! ---------------------------
      ! test
      T_sfc = 0._wp
      
      ! ---------------------------
      ! Go
      ! IF(h_snow.GE.h_Snow_min) THEN   
      !   T_sfc = T_snow                   ! Snow exists, use the snow temperature
      ! ELSE IF(h_ice.GE.h_Ice_min) THEN
      !   T_sfc = T_ice                    ! Ice exists but there is no snow, use the ice temperature
      ! ELSE 
      !   T_sfc = T_wML                    ! No ice-snow cover, use the mixed-layer temperature
      ! END IF
    END SUBROUTINE calc_surface_temperature

#endif
END MODULE mo_lentic_process