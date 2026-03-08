!>
!! @brief Contains parameters for the flake processes
!!
!! @author
!!  Dmitrii Mironov, DWD
!!  Zoe Rehder, MPI-M
!!
!! @par Revision History
!! First version by Dmitrii Mironov (2005/11/17)
!! Adapted for JSBACH by Zoé Rehder (2023/10/10)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_flake_parameters
#ifndef __NO_JSBACH__
  
    USE mo_kind,                ONLY: wp
    USE mo_jsb_impl_constants,  ONLY: def_parameters
  
    IMPLICIT NONE
    PUBLIC
  
    ! ======================================================================================================= !
    !> 1.0 general FLAKE_ parameters
    !>
    ! ..
  
  
    ! ======================================================================================================= !
    !> 2.0 JSBACH specific parameters
    !>
    

    REAL(wp), SAVE :: &
      & depth_w                   = def_parameters,   & ! constant lake depth                   [m] 

      & frac_optic_water(2)       = def_parameters,   & ! Fractions of total radiation flux     [-] 
      & extincoef_optic_water(2)  = def_parameters,   & ! Extinction coefficients per fraction  [m^{-1}]
      & frac_optic_ice(1)         = def_parameters,   & ! Fractions of total radiation flux     [-] 
      & extincoef_optic_ice(1)    = def_parameters,   & ! Extinction coefficients per fraction  [m^{-1}]
      & frac_optic_snow(1)        = def_parameters,   & ! Fractions of total radiation flux     [-] 
      & extincoef_optic_snow(1)   = def_parameters,   & ! Extinction coefficients per fraction  [m^{-1}]

      & albedo_water_ref    = def_parameters,   & ! Water
      & albedo_whiteice_ref = def_parameters,   & ! White ice
      & albedo_blueice_ref  = def_parameters,   & ! Blue ice
      & albedo_bot          = def_parameters,   & ! Albedo of sediment at waterbody bottom [-]

      & C_T_min = def_parameters,   & ! Minimum value of the shape factor C_T (thermocline)
      & C_T_max = def_parameters,   & ! Maximum value of the shape factor C_T (thermocline)
      & C_I_lin = def_parameters,   & ! Shape factor (linear temperature profile in the ice layer)

      & H_Ice_max   = def_parameters,   & ! Maximum ice tickness                            [m] 
      & c_small     = def_parameters,   & ! A small number                                  [-]
      & h_ML_max    = def_parameters,   & ! Maximum mixed-layer depth                       [m]
      & h_Snow_min  = def_parameters,   & ! Minimum snow thickness                          [m]
      & h_Ice_min   = def_parameters,   & ! Minimum ice thickness                           [m]
      & h_ML_min    = def_parameters,   & ! Minimum mixed-layer depth                       [m]
      & u_star_min  = def_parameters,   & ! Minimum value of the surface friction velocity  [m s^{-1}]

      & c_albice_MR = def_parameters,   & ! Constant in the interpolation formula for the ice albedo

      & u_star_min_sf = def_parameters, & ! Minimum value of friction velocity [m s^{-1}]
      & u_wind_min_sf = def_parameters, & ! Minimum wind speed                 [m s^{-1}]
      & c_accur_sf    = def_parameters    ! A small number (accuracy)          [-]
  
    ! ======================================================================================================= !
    !> 3.0 QUINCY specific parameters
    !>
    ! ..

    
    CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_flake_parameters'
  
  CONTAINS
  
    ! the init_constants routines are called in mo_jsb_base:jsbach_setup_models() or :quincy_proc_parameter_init()
  
    ! ======================================================================================================= !
    !> initialize JSBACH specific parameters for the FLALE_ process
    !>
    SUBROUTINE init_flake_parameters
      CHARACTER(len=*), PARAMETER :: routine = modname//':init_flake_parameters'

      ! As a first step, keep lake depths constant
      depth_w   = 10.0_wp  

      ! ! Optical parameters
      ! These parameters should be carefully reviewed before using them
      frac_optic_water(:)      = (/0.10_wp, 0.90_wp/)
      extincoef_optic_water(:) = (/2.0_wp,  0.20_wp/)
      frac_optic_ice(:)        = (/1._wp/)           
      extincoef_optic_ice(:)   = (/1.0E+07_wp/)      
      frac_optic_snow(:)       = (/1._wp/)           
      extincoef_optic_snow(:)  = (/1.0E+07_wp/)      

      albedo_water_ref       = 0.20_wp  
      albedo_whiteice_ref    = 0.60_wp  
      albedo_blueice_ref     = 0.10_wp  
      albedo_bot             = 0.2_wp   
      ! albedo_drysnow_ref     = 0.60_wp  ! Dry snow 
      ! albedo_meltingsnow_ref = 0.10_wp  ! Melting snow 

      ! ! Shape of the thermocline/ice layer
      C_T_min = 0.5_wp  
      C_T_max = 0.9_wp  
      C_I_lin = 0.5_wp  

      ! !  Security constants
      H_Ice_max     = 3._wp     ! Mironov and Ritter (2004MR2004) ice model 
      c_small       = 1.0E-10_wp
      h_ML_max      = 1.0E+3_wp 
      h_Snow_min    = 1.0E-5_wp 
      h_Ice_min     = 1.0E-9_wp
      h_ML_min      = 1.0E-2_wp
      u_star_min    = 1.0E-6_wp
      u_star_min_sf = 1.0E-04_wp
      u_wind_min_sf = 1.0E-02_wp
      c_accur_sf    = 1.0E-07_wp

      !  Empirical parameters.
      c_albice_MR = 95.6_wp ! (Mironov and Ritter 2004)


    END SUBROUTINE init_flake_parameters
    
#endif
END MODULE mo_flake_parameters
  