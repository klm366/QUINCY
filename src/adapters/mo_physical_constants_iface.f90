!>
!! @brief Contains adapters to physcial constants
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-03-18)
!!
MODULE mo_physical_constants_iface
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_physical_constants, ONLY: &
    grav,            & !< Avg. gravitational acceleration        [m/s2]
    stbo,            & !< Stephan-Boltzmann constant             [W/(m2 K4)]
    zemiss_def,      & !< Long-wave surface emissivity factor    []
    rd,              & !< Gas constant                           [J/(K kg)]
    cpd,             & !< Specific heat at constant pressure     [J/(K kg)]
    cvd,             & !< Specific heat at constant volume       [J/(K kg)]
    rhoh2o,          & !< Density of liquid water                [kg/m3]
    alv,             & !< Latent heat for vaporisation           [J/kg]
    als,             & !< Latent heat for sublimation            [J/kg]
    alf,             & !< Latent heat for fusion                 [J/kg]
    tmelt,           & !< Melting temperature of ice/snow        [K]
#ifdef __ICON__
    clw,             & !< Specific heat of liquid water          [J/(K kg)]
    rhoi,            & !< Density of (sea) ice                   [kg/m3]
    rhos,            & !< Density of snow                        [kg/m3]
    ki,              & !< Heat conductivity of ice               [J/(m s K)]
    ks,              & !< Heat conductivity of snow              [J/(m s K)]
#else
    clw => cpliq,    & !< Specific heat of liquid water          [J/(K kg)]
    rhoi => rhoice   & !< Density of (sea) ice                   [kg/m3]
    rhos => rhosno   & !< Density of snow                        [kg/m3]
    ki => alice,     & !< Heat conductivity of ice               [J/(m s K)]
    ks => alsno,     & !< Heat conductivity of snow              [J/(m s K)]
#endif
    ! Auxiliary constants
    rvd1  => vtmpc1, & !< = rv/rd-1                              []
    cpvd1 => vtmpc2, & !< = cpv/cpd-1                            []
    clw                !< Specific heat capacity of liquid water [J/(K kg)]
  IMPLICIT NONE
  PUBLIC

#endif
END MODULE mo_physical_constants_iface
