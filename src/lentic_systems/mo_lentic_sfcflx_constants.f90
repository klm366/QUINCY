!>
!! @brief Contains constants for the processes to compute the surface fluxes from waterbodies
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

MODULE mo_lentic_SfcFlx_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
    
  IMPLICIT NONE
  PUBLIC

     REAL(wp), PARAMETER ::  &
      c_Karman      = 0.40_wp      , & ! The von Karman constant 
      kappa_q_a     = 2.40E-05_wp  , & ! Molecular diffusivity of air for water vapour [m^{2} s^{-1}]
      c_free_conv   = 0.14_wp      , & ! Constant in the expressions for fluxes in free convection
      !  Dimensionless constants in the Monin-Obukhov surface-layer 
      !  similarity relations and in the expressions for the roughness lengths.
      c_MO_u_stab   = 5.0_wp       , & ! Constant of the MO theory (wind, stable stratification)
      c_MO_t_stab   = 5.0_wp       , & ! Constant of the MO theory (temperature, stable stratification)
      c_MO_q_stab   = 5.0_wp       , & ! Constant of the MO theory (humidity, stable stratification)
      c_MO_u_conv   = 15.0_wp      , & ! Constant of the MO theory (wind, convection)
      c_MO_t_conv   = 15.0_wp      , & ! Constant of the MO theory (temperature, convection)
      c_MO_q_conv   = 15.0_wp      , & ! Constant of the MO theory (humidity, convection)
      c_MO_u_exp    = 0.25_wp      , & ! Constant of the MO theory (wind, exponent)
      c_MO_t_exp    = 0.5_wp       , & ! Constant of the MO theory (temperature, exponent)
      c_MO_q_exp    = 0.5_wp       , & ! Constant of the MO theory (humidity, exponent)
      Pr_neutral     = 1.0_wp      , & ! Turbulent Prandtl number at neutral static stability
      Sc_neutral     = 1.0_wp      , & ! Turbulent Schmidt number at neutral static stability
      c_z0u_smooth  = 0.1_wp       , & ! Constant in the expression for z0u (smooth flow) 
      c_z0u_rough    = 1.23E-02_wp , & ! The Charnock constant in the expression for z0u (rough flow)
      c_z0u_rough_L = 1.00E-01_wp  , & ! An increased Charnock constant (used as the upper limit)
      c_z0u_ftch_f  = 0.70_wp      , & ! Factor in the expression for fetch-dependent Charnock parameter
      c_z0u_ftch_ex = 0.3333333_wp , & ! Exponent in the expression for fetch-dependent Charnock parameter
      c_z0t_rough_1 = 4.0_wp       , & ! Constant in the expression for z0t (factor) 
      c_z0t_rough_2 = 3.2_wp       , & ! Constant in the expression for z0t (factor)
      c_z0t_rough_3 = 0.5_wp       , & ! Constant in the expression for z0t (exponent) 
      c_z0q_rough_1 = 4.0_wp       , & ! Constant in the expression for z0q (factor)
      c_z0q_rough_2 = 4.2_wp       , & ! Constant in the expression for z0q (factor)
      c_z0q_rough_3 = 0.5_wp       , & ! Constant in the expression for z0q (exponent)
      c_z0t_ice_b0s = 1.250_wp     , & ! Constant in the expression for z0t over ice
      c_z0t_ice_b0t = 0.149_wp     , & ! Constant in the expression for z0t over ice
      c_z0t_ice_b1t = -0.550_wp    , & ! Constant in the expression for z0t over ice
      c_z0t_ice_b0r = 0.317_wp     , & ! Constant in the expression for z0t over ice
      c_z0t_ice_b1r = -0.565_wp    , & ! Constant in the expression for z0t over ice
      c_z0t_ice_b2r = -0.183_wp    , & ! Constant in the expression for z0t over ice
      c_z0q_ice_b0s = 1.610_wp     , & ! Constant in the expression for z0q over ice
      c_z0q_ice_b0t = 0.351_wp     , & ! Constant in the expression for z0q over ice
      c_z0q_ice_b1t = -0.628_wp    , & ! Constant in the expression for z0q over ice
      c_z0q_ice_b0r = 0.396_wp     , & ! Constant in the expression for z0q over ice
      c_z0q_ice_b1r = -0.512_wp    , & ! Constant in the expression for z0q over ice
      c_z0q_ice_b2r = -0.180_wp    , & ! Constant in the expression for z0q over ice
      Re_z0s_ice_t  = 2.5_wp       , & ! Threshold value of the surface Reynolds number 
                                        ! used to compute z0t and z0q over ice (Andreas 2002)
      ! Re_z0u_thresh = 0.1_wp          ! Threshold value of the roughness Reynolds number 
                                        ! [value from Zilitinkevich, Grachev, and Fairall (200),
                                        ! currently not used] 
      z0u_ice_rough = 1.0E-03_wp        ! Aerodynamic roughness of the ice surface [m] (rough flow)
      
    CHARACTER(len=*), PARAMETER, PRIVATE :: modname =  'mo_lentic_SfcFlx_constants'
  
#endif     
END MODULE mo_lentic_SfcFlx_constants