!>
!! @brief Contains constants for the radiation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-15)
!!
MODULE mo_radiation_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::                  &
     Epar                   = 2.2E5_wp,   & ! Energy content of PAR [J / mol(photons)]=(4.6 mol/MJ PAR)**-1.  R: for update_par                                     
     SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region.  R: for update_par  
!     LaiMax = 8._wp,                      & ! Maximum LAI (used for nitrogen scaling) 
     LaiMin = 1.E-9_wp,                   & ! Minimum Lai in PAR computation.R: for calc_fapar_lai_cl
     LaiLimit = 3._wp,                    & ! Minimum LAI: used for nitrogen scaling and estimation of fractional cover. R: for calc_fapar_lai_cl
     FcMax = 1.0_wp,                      & ! Maximum fractional vegetation cover. R: for calc_fapar_lai_cl
     FcMin = 1.E-3_wp,                    & ! Minimum fractional vegetation cover. R: for calc_fapar_lai_cl
     ZenithMinPar = 1.E-3_wp                ! Minimum cos of zenith angle for which PAR is calculated. R: for calc_fapar_lai_cl
!     ZenithMin = 0.0174524_wp,            & ! Check for solar zenith angle > 89 degrees
!     SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region
!  
!     & minStomaConductance = 0.0_wp,        & ! Minimum stomatal conductance [mol H2O /(m^2 s) ??]

#endif
END MODULE mo_radiation_constants
