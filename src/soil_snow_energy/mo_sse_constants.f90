!>
!! @brief Contains constants for the soil and snow energy process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_sse_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_physical_constants,  ONLY: ks, ki, clw, cs, ci, rhoh2o, dens_snow, dens_ice

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::                    &
    ! Parameters for solid soil particulate
    & hcond_dry_soil     = 0.05_wp,         & !< Dry soil heat conductivity (O'Donnel et al 2009)             [W/m/K]
    ! Parameters for liquid water
    & vol_hcap_water     = clw * rhoh2o,    & !< Volumetric heat capacity of liquid water                     [J/m3/K]
    & hcond_water        = 0.57_wp,         & !< Thermal conductivity of liquid water (Hillel 1982)           [W/K/m]
    ! Parameters for frozen water in soil (ice)
    ! @todo: vol_hcap_ice is 2090000 in jsbach3
    & vol_hcap_ice       = ci * dens_ice,   & !< =1931202 Volumetric heat capacity of ice                     [J/m3/K]
    & hcond_ice          = ki,              & !< =2.1656  thermal conductivity of ice                         [W/K/m]
    ! Parameters for organic soil component
    & vol_hcap_org       = 2.5e+06_wp,      & !< Organic layer volumetric heat capacity (Beringer et al 2001) [J/m3/K]
    & hcond_org          = 0.25_wp,         & !< Organic layer heat conductivity (Beringer et al 2001)        [W/m/K]
    & hcond_dry_org      = 0.05_wp,         & !< Dry soil heat conductivity (Farouki 1981)                    [W/m/K]
    ! Parameters for bedrock
    & vol_hcap_bedrock   = 2.e+06_wp,       & !< Volumetric heat capacity of bedrock (Bonan 2002)             [J/m3/K]
    & hcond_bedrock      = 2._wp,           & !< Thermal conductivity of bedrock (Tarnocai 2004)              [W/m/K]
    ! Parameters for free soil pores (air)
    & vol_hcap_air       = 1300._wp,        & !< Volumetric heat capacity of free soil pores (from VIC model) [J/m3/K]
    ! Parameters for snow
    & snow_depth_min     = 0.02_wp,         & !< Minimum snow depth below which no snow layer is considered   [m]
    ! @todo: vol_hcap_snow is 634500 in jsbach3
    & vol_hcap_snow      = cs * dens_snow,  & !< =627000 Volumetric heat capacity of snow                     [J/m3/K]
    & vol_hcap_snow_min  = 526500._wp,      & !< Minimum snow heat capacity (fresh snow)                      [J/m3/K]
    & hcond_snow         = ks,              & !< =0.31 Thermal conductivity of snow                           [W/m/K]
    ! @todo: hcond_snow_min is 0.1 in jsbach3
    & hcond_snow_min     = 0.03_wp            !< Minimum thermal conductivity fo snow (fresh snow)            [W/m/K]

  !! quincy
  REAL(wp), SAVE :: &
    soil_heat_cap,  &           !< soil heat capacity (temporarily for SLM offline) J/kg/K 
    !rho_bulk_soil,  &           !< soil bulk denisty (kg / m3) 
    soil_therm_cond, &          !< soil thermal conductivity (W m-1 K-1)
    ! constants in simple multi-layer soil moisture scheme to describe the 
    InterceptionEfficiency, &   !< Efficiency of interception of precipitation as rain
    w_skin_max, &               !< maximum water storage per unit LAI (m)
    drain_frac, &               !< preferential flow fraction of infiltrating water per meter soil 
    fdrain_runoff               !< preferential flow fraction of runoff

  REAL(wp), SAVE :: &
    height_wind                 !< ...

 ! constants in pedo-transfer functions according to Saxton & Rawls 2006
 REAL(wp) :: & 
    wpot_pwp, &                 !< water potential at permanent wilitng point (Pa)
    wpot_fc, &                  !< water potential at field capacity (Pa)
    wpot_min, &                 !< minimum water potential (MPa) for numerical reasons
    k_pwp_a, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_c, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_s, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_sc, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_at, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    k_pwp_bt, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    k_fc_a, &                   !< empirical constant in pedo-transfer function to get water content at FC 
    k_fc_c, &                   !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_s, &                   !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_sc, &                  !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_at, &                  !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_bt, &                  !< empirical constant in pedo-transfer function to get water content at FC
    k_fc_ct, &                  !< empirical constant in pedo-transfer function to get water content at FC
    k_sat_a, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_c, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_s, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_sc, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_at, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_bt, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_ct, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    k_sat_dt, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    kdiff_sat_max               !< maximum diffusivity of soil at saturation (m/s)

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sse_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> quincy: initialize constants from sse 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_sse_constants

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_sse_constants'

    soil_heat_cap           = 1480.0_wp     ! textbook for wet sandy soil
    !rho_bulk_soil           = 1695.86_wp    ! Bernhard Ahrens for Hainich -> using site specific value now
    soil_therm_cond         = 2.0_wp        ! Text book for wet sandy soil
    InterceptionEfficiency  = 0.25_wp       ! JSBACH4 
    w_skin_max              = 0.0002_wp     ! JSBACH4 
    drain_frac              = 0.01_wp       ! tuned
    fdrain_runoff           = 0.95_wp       ! set from ORCHIDEE
    height_wind             = 10.0_wp       ! JSBACH ?

    ! Saxton & Rawls, 2006, includes conversion from fractional to mass% texture
    wpot_pwp        = 1500._wp
    wpot_fc         = 33.0_wp
    wpot_min        = -15.0_wp  ! no plant can survive 8 MPa...
    k_pwp_a         = 0.031_wp 
    k_pwp_c         = 0.487_wp 
    k_pwp_s         = -0.024_wp 
    k_pwp_sc        = 0.068_wp 
    k_pwp_at        = -0.02_wp
    k_pwp_bt        = 0.14_wp
    k_fc_a          = 0.299_wp
    k_fc_c          = 0.195_wp 
    k_fc_s          = -0.251_wp 
    k_fc_sc         = 0.452_wp 
    k_fc_at         = -0.015_wp
    k_fc_bt         = -0.373_wp
    k_fc_ct         = 1.293_wp
    k_sat_a         = 0.078_wp
    k_sat_c         = 0.034_wp 
    k_sat_s         = 0.278_wp 
    k_sat_sc        = -0.584_wp 
    k_sat_at        = -0.107_wp
    k_sat_bt        = 0.6360_wp
    k_sat_ct        = 0.097_wp 
    k_sat_dt        = 0.043_wp
    kdiff_sat_max   = 1930.0_wp / 1000.0_wp / 3600._wp ! m/s

  END SUBROUTINE init_sse_constants

#endif
END MODULE mo_sse_constants
