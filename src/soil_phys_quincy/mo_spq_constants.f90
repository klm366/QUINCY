!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains constants for the spq process
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-21)
!!
MODULE mo_spq_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC


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
    height_wind                 !< 

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


  ! constants used as factor for multiplication with namelist parameters (spq_config) when "parameter sensitivity" is enabled
  ! see spq_config & mo_quincy_set_parameters:Modify_quincy_namelist_parameter
  REAL(wp) :: & 
    f_psensi_soil_sand, &
    f_psensi_soil_silt, &
    f_psensi_soil_awc_prescribe, &
    f_psensi_bulk_density


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_spq_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from spq 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_spq_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    ! 
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_spq_constants'


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


    ! constants used as factor for multiplication with namelist parameters (spq_config) when "parameter sensitivity" is enabled
    f_psensi_soil_sand            = 1.0_wp
    f_psensi_soil_silt            = 1.0_wp
    f_psensi_soil_awc_prescribe   = 1.0_wp
    f_psensi_bulk_density         = 1.0_wp

  END SUBROUTINE init_spq_constants

#endif
END MODULE mo_spq_constants
