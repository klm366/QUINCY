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
    argas,           & !< Molar/universal/ideal gas constant     [J/(K mol)]
    zemiss_def,      & !< Long-wave surface emissivity factor    []
    rd,              & !< Gas constant                           [J/(K kg)]
    cpd,             & !< Specific heat at constant pressure     [J/(K kg)]
    cvd,             & !< Specific heat at constant volume       [J/(K kg)]
    rhoh2o,          & !< Density of liquid water                [kg/m3]
    alv,             & !< Latent heat for vaporisation           [J/kg]
    als,             & !< Latent heat for sublimation            [J/kg]
    alf,             & !< Latent heat for fusion                 [J/kg]
    tmelt,           & !< Melting temperature of ice/snow        [K]
    amd,             & !< Molar weight of dry air                [g/mol]
    amco2,           & !< Molar weight of CO2                    [g/mol]
#ifdef __ICON__
    clw,             & !< Specific heat of liquid water          [J/(K kg)]
    ci,              & !< Heat capacity of ice                   [J/(K kg)]
    cs,              & !< Specific heat of snow                  [J/(K kg)]
    rhoi,            & !< Density of (sea) ice                   [kg/m3]
    rhos,            & !< Density of snow                        [kg/m3]
    ki,              & !< Heat conductivity of ice               [J/(m s K)]
    ks,              & !< Heat conductivity of snow              [J/(m s K)]
#else
    clw => cpliq,    & !< Specific heat of liquid water          [J/(K kg)]
    ci => cpice,     & !< Specific heat of ice                   [J/(K kg)]
    cs => cpsno,     & !< Specific heat of snow                  [J/(K kg)]
    rhoi => rhoice,  & !< Density of (sea) ice                   [kg/m3]
    rhos => rhosno,  & !< Density of snow                        [kg/m3]
    ki => alice,     & !< Heat conductivity of ice               [J/(m s K)]
    ks => alsno,     & !< Heat conductivity of snow              [J/(m s K)]
#endif
    ! Auxiliary constants
    rvd1  => vtmpc1, & !< = rv/rd-1                              []
    cpvd1 => vtmpc2    !< = cpv/cpd-1                            []

! cvdifts in echam is not a PARAMETER, therefore it can't be USEd from modules here
! TBD: Use cvdifts=1 for standalone model?
!#ifndef STANDALONE
!#ifdef __ICON__
!  USE mo_echam_vdiff_params, ONLY: cvdifts
!#else
!  USE mo_physc2, ONLY: cvdifts
!#endif
!#endif

  IMPLICIT NONE
  PUBLIC

!#ifdef STANDALONE
!  REAL(wp), PARAMETER :: cvdifts = 1._wp
!#else
  REAL(wp), PARAMETER :: cvdifts = 1.5_wp
!#endif
  REAL(wp), PARAMETER :: tpfac1 = cvdifts
  REAL(wp), PARAMETER :: tpfac2 = 1._wp / tpfac1
  REAL(wp), PARAMETER :: tpfac3 = 1._wp - tpfac2 

#endif
END MODULE mo_physical_constants_iface
