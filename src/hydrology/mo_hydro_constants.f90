!>
!! @brief Contains constants for the radiation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-15)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hydro_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::                   &
    ! Parameters used for snow cover fraction
    & wsn2fract_eps    = 1.E-12_wp,        &
    & wsn2fract_sigfac = 0.15_wp,          &
    & wsn2fract_const  = 0.95_wp,          & !< Inverse of snow height when snow is considered to completely cover the ground

    & InterceptionEfficiency = 0.25_wp,    & !< Efficiency of interception of precipitation (rain and snow)

    & Epar                   = 2.2E5_wp,   & ! Energy content of PAR [J / mol(photons)]=(4.6 mol/MJ PAR)**-1.  R: for update_par
    & SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region.  R: for update_par
!    & LaiMax = 8._wp,                      & ! Maximum LAI (used for nitrogen scaling)
    & FcMax = 1.0_wp,                      & ! Maximum fractional vegetation cover. R: for calc_fapar_lai_cl
    & FcMin = 1.E-3_wp,                    & ! Minimum fractional vegetation cover. R: for calc_fapar_lai_cl
    & ZenithMinPar = 1.E-3_wp,             &  ! Minimum cos of zenith angle for which PAR is calculated. R: for calc_fapar_lai_cl
!    & ZenithMin = 0.0174524_wp,            & ! Check for solar zenith angle > 89 degrees
!    & SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region
!
!    & minStomaConductance = 0.0_wp,        & ! Minimum stomatal conductance [mol H2O /(m^2 s) ??]

    ! Parameters describing orographic steepness in infiltration computation
    & oro_var_min = 100._wp,               &
    & oro_var_max = 1000._wp,              &

    ! Parameters for the computation of canopy conductance/resistance using Eq. 3.3.2.12 in ECHAM3 manual
    & conductance_k = 0.9_wp,              & !<
    & conductance_a = 5000._wp,            & !< Jm-3
    & conductance_b = 10._wp,              & !< Wm-2
    & conductance_c = 100._wp,             & !< ms-1

    ! Parameters for organic soil component (top layer)
    & vol_porosity_org_top    = 0.95_wp,       & !< volumetric porosity of organic part of soil layer [m/m]
    & vol_field_cap_org_top   = 0.95_wp,       & !< volumetric field capacity of organic layer [m/m]
    & vol_p_wilt_org_top      = 0.255_wp,      & !< volumetricwilting point of organic layer [m/m]
    & hyd_cond_sat_org_top    = 0.0001_wp,     & !< Saturated hydraulic conductivity of organic part of soil layer [m/s]
    & bclapp_org_top          = 4._wp,         & !< Exponent b in Clapp and Hornberger of organic part of soil layer  []
    & matrix_pot_org_top      = 0.1_wp,        & !< Matrix potential of organic part of soil layer [m]
                                                 ! (values (roughy) from Chadburn 2022)
    & pore_size_index_org_top = 0.7_wp,        & !< Pore size distribution index of organic part of soil layer []
    ! Parameters for organic soil component (deep layers)
    & vol_porosity_org_below  = 0.82_wp,       & !< volumetric porosity of organic part of soil layer [m/m]
    & vol_field_cap_org_below = 0.82_wp,       & !< volumetric field capacity of organic layer [m/m]
    & vol_p_wilt_org_below    = 0.255_wp,      & !< volumetricwilting point of organic layer [m/m]
    & hyd_cond_sat_org_below  = 0.00000001_wp, & !< Saturated hydraulic conductivity of organic part of soil layer [m/s]
    & bclapp_org_below        = 8._wp,         & !< Exponent b in Clapp and Hornberger of organic part of soil layer  []
    & matrix_pot_org_below    = 1.0_wp,        & !< Matrix potential of organic part of soil layer [m]
                                                 ! (values (roughy) from Chadburn 2022)
    & pore_size_index_org_below = 0.7_wp,      & !< Pore size distribution index of organic part of soil layer []
    & thresh_org          = 0.5_wp,            & !< Threshold of organic material above which connected flow pathways form

    & beta_perc           = 0.139_wp             !< Parameter from percolation theory

  ! quincy_jsbach: will be removed and re-implemented >>>
  ! constants in pedo-transfer functions according to Saxton & Rawls 2006
  REAL(wp) :: & 
    & wpot_pwp, &                 !< water potential at permanent wilitng point (Pa)
    & wpot_fc, &                  !< water potential at field capacity (Pa)
    & wpot_min, &                 !< minimum water potential (MPa) for numerical reasons
    & k_pwp_a, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    & k_pwp_c, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    & k_pwp_s, &                  !< empirical constant in pedo-transfer function to get water content at PWP
    & k_pwp_sc, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    & k_pwp_at, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    & k_pwp_bt, &                 !< empirical constant in pedo-transfer function to get water content at PWP
    & k_fc_a, &                   !< empirical constant in pedo-transfer function to get water content at FC 
    & k_fc_c, &                   !< empirical constant in pedo-transfer function to get water content at FC
    & k_fc_s, &                   !< empirical constant in pedo-transfer function to get water content at FC
    & k_fc_sc, &                  !< empirical constant in pedo-transfer function to get water content at FC
    & k_fc_at, &                  !< empirical constant in pedo-transfer function to get water content at FC
    & k_fc_bt, &                  !< empirical constant in pedo-transfer function to get water content at FC
    & k_fc_ct, &                  !< empirical constant in pedo-transfer function to get water content at FC
    & k_sat_a, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    & k_sat_c, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    & k_sat_s, &                  !< empirical constant in pedo-transfer function to get water content at SAT
    & k_sat_sc, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    & k_sat_at, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    & k_sat_bt, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    & k_sat_ct, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    & k_sat_dt, &                 !< empirical constant in pedo-transfer function to get water content at SAT
    & kdiff_sat_max               !< maximum diffusivity of soil at saturation (m/s)


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sse_constants'

CONTAINS

  ! ======================================================================================================= !
  !>quincy: initialize constants from HYDRO_ - originally from SSE_
  !> 
  SUBROUTINE init_hydro_constants

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_hydro_constants'

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

  END SUBROUTINE init_hydro_constants
  ! quincy_jsbach<<<

#endif
END MODULE mo_hydro_constants
