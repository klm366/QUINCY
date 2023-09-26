!>
!! @brief Contains constants for the lentic processes
!!
!! @author
!!  Dmitrii Mironov, DWD
!!  Zoe Rehder, MPI-M
!!
!! @par Revision History
!! First version by Dmitrii Mironov (2005/11/17)
!! Adapted for JSBACH by Zoé Rehder (2023/02/09)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_lentic_constants
#ifndef __NO_JSBACH__

USE mo_kind, ONLY: wp

IMPLICIT NONE
PUBLIC

     REAL(wp), PARAMETER ::  &
      & depth_w   = 10.0_wp       , & ! constant lake depth [m] ! As a first step, keep lake depths constant
      & rot_earth = 7.2921E-05_wp , & ! rotation rate of the earth  [rad/s]


      !  Dimensionless constants 
      !  in the equations for the mixed-layer depth 
      !  and for the shape factor with respect to the temperature profile in the thermocline
      & c_cbl_1       = 0.17_wp            , & ! Constant in the CBL entrainment equation
      & c_cbl_2       = 1._wp              , & ! Constant in the CBL entrainment equation
      & c_sbl_ZM_n    = 0.5_wp             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
      & c_sbl_ZM_s    = 10._wp             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
      & c_sbl_ZM_i    = 20._wp             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
      & c_relax_h     = 0.030_wp           , & ! Constant in the relaxation equation for the SBL depth
      & c_relax_C     = 0.0030_wp          , & ! Constant in the relaxation equation for the shape factor
                                                  ! with respect to the temperature profile in the thermocline

      & albedo_bot    = 0.2_wp             , & ! Albedo of sediment at waterbody bottom [-]
      & C_T_min       = 0.5_wp             , & ! Minimum value of the shape factor C_T (thermocline)
      & C_T_max       = 0.9_wp             , & ! Maximum value of the shape factor C_T (thermocline)
      & C_I_lin       = 0.5_wp             , & ! Shape factor (linear temperature profile in the ice layer)
      & C_I_MR        = 1._wp/12._wp       , & ! Constant in the MR2004 expression for I shape factor

      & Phi_I_pr0_lin = 1._wp              , & ! I shape-function derivative (linear profile) 
      & Phi_I_pr1_lin = 1._wp              , & ! I shape-function derivative (linear profile) 
      & Phi_I_ast_MR  = 2._wp              , & ! Constant in the MR2004 expression for I shape factor
      & H_Ice_max     = 3._wp              , & ! Maximum ice tickness in 
                                               ! the Mironov and Ritter (2004, MR2004) ice model [m] 
      & C_TT_1        = 11._wp/18._wp      , & ! Constant in the expression for C_TT (thermocline)
      & C_TT_2        = 7._wp/45._wp       , & ! Constant in the expression for C_TT (thermocline)

      ! ! Thermodynamic parameters
      ! & grav = tpl_grav          = 9.81_wp       , & ! Acceleration due to gravity [m s^{-2}]
      & tpl_T_r           = 277.13_wp     , & ! Temperature of maximum density of fresh water [K]
      & tpl_a_T           = 1.6509E-05_wp , & ! Constant in the fresh-water equation of state [K^{-2}]

      & tpl_c_I           = 2.1E+03_wp    , & ! Specific heat of ice [J kg^{-1} K^{-1}]

      & tpl_kappa_I       = 2.29_wp       , & ! Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
      & tpl_kappa_w       = 5.46E-01_wp   , & ! Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
      !! & alf = tpl_L_f           = 3.3E+05_wp    , & ! Latent heat of fusion [J kg^{-1}]

      & tpl_T_f           = 273.15_wp     , & ! Fresh water freezing point [K]
      ! & rhoice = tpl_rho_I         = 9.1E+02_wp    , & ! Density of ice [kg m^{-3}]
      & rhoice            =  916.7_wp     , & ! density of pure ice [kg/m3]  , !! copied from ICON
      & con_m             = 1.50E-5_wp    , & ! [m^2/s]  kinematic viscosity of dry air, !! copied from ICON
      & con_h             = 2.20E-5_wp    , & ! [m^2/s]  scalar conductivity of dry air, !! copied from ICON
      & tpl_rho_w_r       = 1.0E+03_wp    , & ! Maximum density of fresh water [kg m^{-3}]

      ! !  Parameters of the shape functions 
      ! !  Indices refer to T - thermocline, S - snow, I - ice,
      ! !  B1 - upper layer of the bottom sediments, B2 - lower layer of the bottom sediments.
      ! !  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
      ! !  at "zeta=0" ad "zeta=1", respectively.
      & Phi_T_pr0_1   = 40._wp/3._wp   , & ! Constant in the expression for the T shape-function derivative 
      & Phi_T_pr0_2   = 20._wp/3._wp   , & ! Constant in the expression for the T shape-function derivative 

      ! !  Security constants
      & c_small        = 1.0E-10_wp, & ! A small number
      & h_ML_max       = 1.0E+3_wp , & ! Maximum mixed-layer depth [m]
      & h_Snow_min     = 1.0E-5_wp , & ! Minimum snow thickness [m]
      & h_Ice_min      = 1.0E-9_wp , & ! Minimum ice thickness [m]
      & h_ML_min       = 1.0E-2_wp , & ! Minimum mixed-layer depth [m]
      & u_star_min     = 1.0E-6_wp , & ! Minimum value of the surface friction velocity [m s^{-1}]

      ! !  Albedo for water, ice and snow.
      & albedo_water_ref       = 0.20_wp  , & ! Water
      & albedo_whiteice_ref    = 0.60_wp  , & ! White ice
      & albedo_blueice_ref     = 0.10_wp  , & ! Blue ice
      ! & albedo_drysnow_ref     = 0.60_wp  , & ! Dry snow 
      ! & albedo_meltingsnow_ref = 0.10_wp  , & ! Melting snow 

      !  Empirical parameters.
      & c_albice_MR = 95.6_wp, &        ! Constant in the interpolation formula for 
                                      ! the ice albedo (Mironov and Ritter 2004)

      ! tpl_c_w = clw
      ! These parameters should be carefully reviewed before using them
      & frac_optic_water(2)      = (/0.10_wp, 0.90_wp/),  &! Fractions of total radiation flux [-] 
      & extincoef_optic_water(2) = (/2.0_wp,  0.20_wp/),  &! Extinction coefficients per fraction [m^{-1}]
      & frac_optic_ice(1)        = (/1._wp/),             &! Fractions of total radiation flux [-] 
      & extincoef_optic_ice(1)   = (/1.0E+07_wp/),        &! Extinction coefficients per fraction [m^{-1}]
      & frac_optic_snow (1)      = (/1._wp/),             &! Fractions of total radiation flux  [-] 
      & extincoef_optic_snow(1)  = (/1.0E+07_wp/)          ! Extinction coefficients per fraction [m^{-1}]
        
      CHARACTER(len=*), PARAMETER, PRIVATE :: modname =  'mo_lentic_constants'
#endif 
  
END MODULE mo_lentic_constants
