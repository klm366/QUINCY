!> Contains constants used in the hydrology processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Constants and parameters used in the hydrology process
!>
!> The module defines physical parameters needed in hydrological parametrisations as well as more
!> technical parameters e.g. to ease the handling of namelist parameters.
!>
!> _References_
!>
!> Parametrization of snow fraction
!>
!> - A. Roesch et al. Climate Dynamics (2001) 17: 933-947
!>
!> Parametrization of runoff/drainage
!>
!> - Todini, E., The ARNO rainfall—runoff model, Journal of Hydrology 175 (1996) 339-382
!> - Dümenil, L., & Todini, E. (1992). A rainfall-runoff scheme for use in the Hamburg climate model.
!>   In J. P. O'Kane (Ed.), Advances in Theoretical Hydrology: A tribute to James Dooge (pp. 129-157).
!>   Amsterdam: Elsevier Science Publishers B.V.
!>
!> Parametrization of canopy conductance/resistance
!>
!> - The ECHAM 3 Atmospheric General Circulation Model, technical report No.6,
!>   Deutsches Klimarechnezentrum (1994)
!>
!> Parametrization of organic soil components
!>
!>  - Letts et al., 2000
!>  - Chadburn 2022
!>  - Clapp and Hornberger
!>
!> Hydaulic conductivity
!>
!> - Freeze, R. A. and Cherry, J. A. (1979): Groundwater, ISBN: 0-13-365312-9
!>
!>
MODULE mo_hydro_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER :: &
    ! Identifiers used for the pedotransfer functions of a specific soil hydrology model
    & VanGenuchten_     = 1,                & !< Identifier for VanGenuchten soil hydrology model
    & BrooksCorey_      = 2,                & !< Identifier for Brooks and Corey soil hydrology model
    & Campbell_         = 3,                & !< Identifier for Campbell soil hydrology model

    ! Identifiers used for the interpolation scheme
    & Upstream_         = 1,                & !< Identifier to use K from upper soil layer and MAX(D) of upper or actual layer
    & Arithmetic_       = 2,                & !< Identifier to average K and D between actual and upper soil layer

    ! Identifiers used for the pond dynamics scheme
    & Quad_             = 1,                & !< Identifier for quadratic scaling function between pond depth and area
    & Tanh_             = 2,                & !< Identifier for tanh scaling function between pond depth and area

    ! Identifiers used for hydrological scale parametrizations
    & Semi_Distributed_ = 1,                & !< Identifier to use ARNO scheme for soil hydrology
    & Uniform_          = 2                   !< Identifier to assume uniform grid cell characteristics

  REAL(wp), PARAMETER :: &
    ! Parameters used to calculate snow fraction from snow amount
    & wsn2fract_eps    = 1.E-12_wp,        & !< Used for numerical reasons (Roesch et al. 2001, Eq. 6)
    & wsn2fract_sigfac = 0.15_wp,          & !< Factor to approximate sub-grid scale slopes (Roesch et al. 2001, Eq. 7)
    & wsn2fract_const  = 0.95_wp,          & !< Factor for snow fraction (Roesch et al. 2001, Eq. 6 and 7)

    ! Canopy interception
    & InterceptionEfficiency = 0.25_wp,    & !< Efficiency of precipitation interception (rain and snow)

    & Epar                   = 2.2E5_wp,   & !TODO: remove, not used in HYDRO
    & SoilReflectivityParMin = 0.0_wp,     & !TODO: remove, not used in HYDRO

    & FcMax = 1.0_wp,                      & !TODO: remove, not used in HYDRO
    & FcMin = 1.E-3_wp,                    & !TODO: remove, not used in HYDRO
    & ZenithMinPar = 1.E-3_wp,             & !TODO: remove, not used in HYDRO

    ! Parameters required for ARNO Scheme
    & oro_var_min = 100._wp,               & !< ARNO Scheme minimum orographic standard deviation threshold
    & oro_var_max = 1000._wp,              & !< ARNO Scheme maximum orographic standard deviation threshold
    & drain_min   = 0.001_wp / (3600._wp*1000._wp), & !< maximum flux for slow subsurface drainage, i.e. below field capacity
                                                      !< (Todini, 1996; Duemenil, 1992)
    & drain_max   = 0.1_wp   / (3600._wp*1000._wp), & !< maximum flux for fast subsurface drainage, i.e. above field capacity
                                                      !< (Todini, 1996; Duemenil, 1992)
    & drain_exp   = 1.5_wp,                & !< drainage scaling exponent for ARNO scheme

    ! Parameters used for pond computations
    & oro_crit = 100._wp,                  & !< Reference topographic standard deviation for pond scaling factor [m]

    ! Parameters for the computation of canopy conductance/resistance using Eq. 3.3.2.12 in ECHAM3 manual
    & conductance_k = 0.9_wp,              & !< Parameter for canopy conductance/resistance
    & conductance_a = 5000._wp,            & !< Parameter for conductance/resistance [Jm-3]
    & conductance_b = 10._wp,              & !< Parameter for canopy conductance/resistance [Wm-2]
    & conductance_c = 100._wp,             & !< Parameter for conductance/resistance [ms-1]

    ! Parameters for organic soil component (top layer)
    & vol_porosity_org_top    = 0.95_wp,       & !< Volumetric porosity of (top) organic layer [m/m]
    & vol_field_cap_org_top   = 0.95_wp,       & !< Volumetric field capacity of (top) organic layer [m/m]
    & vol_p_wilt_org_top      = 0.255_wp,      & !< Volumetric wilting point of (top) organic layer [m/m]
    & vol_wres_org_top        = 0.050_wp,      & !< Volumetric residual water content of (top) organic layer [m/m] (Letts et al., 2000)
    & hyd_cond_sat_org_top    = 0.0001_wp,     & !< Saturated hydraulic conductivity of (top) organic part of soil layer [m/s]
    & bclapp_org_top          = 4._wp,         & !< Exponent b in Clapp and Hornberger of (top) organic part of soil layer  []
    & matric_pot_org_top      = -0.1_wp,       & !< Soil matric potential of (top) organic layer [m]
                                                 !< (values (roughy) from Chadburn 2022)
    & pore_size_index_org_top = 0.7_wp,        & !< Pore size distribution index of (top) organic layer []

    ! Parameters for organic soil component (deep layers)
    & vol_porosity_org_below  = 0.82_wp,       & !< Volumetric porosity of (deep) organic layers [m/m]
    & vol_field_cap_org_below = 0.82_wp,       & !< Volumetric field capacity of (deep) organic layers [m/m]
    & vol_p_wilt_org_below    = 0.255_wp,      & !< Volumetric wilting point of (deep) organic layers [m/m]
    & vol_wres_org_below      = 0.150_wp,      & !< Volumetric residual water content of (deep) organic layers [m/m] (Letts et al., 2000)
    & hyd_cond_sat_org_below  = 0.00000001_wp, & !< Saturated hydraulic conductivity of (deep) organic layers [m/s]
    & bclapp_org_below        = 8._wp,         & !< Exponent b in Clapp and Hornberger of (deep) organic layers  []
    & matric_pot_org_below    = -1.0_wp,       & !< Soil matric potential of (deep) organic layers [m]
                                                 !< (values (roughy) from Chadburn 2022)
    & pore_size_index_org_below = 0.7_wp,      & !< Pore size distribution index of (deep) organic layers []
    & thresh_org              = 0.5_wp,        & !< Threshold of organic material above which connected flow pathways form
    & beta_perc               = 0.139_wp,      & !< Parameter from percolation theory

    ! Parameters for bedrock from  Freeze and Cherry, 1979: Groundwater
    & k_brock                 = 1.E-6_wp         !< Hydraulogic conductivity (fractured) bedrock [m/s] (Freeze & Cherry, 1979)

#endif
END MODULE mo_hydro_constants
