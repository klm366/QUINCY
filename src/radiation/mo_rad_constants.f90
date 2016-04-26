!>
!! @brief Contains constants for the radiation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-15)
!!
MODULE mo_rad_constants
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
     ZenithMinPar = 1.E-3_wp,             & ! Minimum cos of zenith angle for which PAR is calculated. R: for calc_fapar_lai_cl
!     ZenithMin = 0.0174524_wp,            & ! Check for solar zenith angle > 89 degrees
!     SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region
!
!     & minStomaConductance = 0.0_wp,        & ! Minimum stomatal conductance [mol H2O /(m^2 s) ??]
    ! Parameters for spatially constant albedo of the soil surface (UseAlbedoSoilConst)
    AlbedoSoilConstVis = 0.2,             & ! Albedo constant alternatively usable for alb_vis_mineralsoil
    AlbedoSoilConstNir = 0.2,             & ! See above, for the nir range
    !
    SkyViewFactor      = 0.5,             & ! Constant in calculating the sky view fractions (only to scale the e-function)
    ZenithAngleFactor  = 2.,              & ! Factor in solar zenith angle dependence of snow albedo
    !                                 ! (the increase of snow albedo is the higher this factor
    ! Parameters of snow age scheme albedo
    AlbedoSnowVisMax   = 0.90,            & ! Maximum albedo of fresh snow in the visible range
    AlbedoSnowNirMax   = 0.60,            & ! Maximum albedo of fresh snow in the NIR range
    AlbedoSnowVisAge   = 0.15,            & ! Maximal rel. reduction of snow albedo by aging in the visible range
    AlbedoSnowNirAge   = 0.5,             & ! Maximal rel. reduction of snow albedo by aging in the NIR range
    AlbedoCanopySnow_age = 0.20,          &  ! Albedo of snow covered canopy for snow age scheme (BATS)
    AlbedoSnowAngle    = 0.4,             & ! Maximal rel. reduction of snow absorption by large solar zenith angle
    !
    ! Parameters of snow temperature scheme albedo
    TempAlbedoSnowMax     = 5.0,          & ! Maximum snow albedo at this temperature below melting point of H2O
    AlbedoCanopySnow_temp = 0.20,         & ! Albedo of snow covered canopy for temperature scheme
    !
    ! Parameters of glacier scheme albedo
    AlbedoGlacierVisMin  = 0.78,          & ! Albedo of glacier in the visible range at the melting point
    AlbedoGlacierVisMax  = 0.9,           & ! Albedo of glacier in the visible range at hard frost
    TempAlbedoGlacierMax = 5.0,           & ! Maximum glacier albedo at this temperature below melting point of H2O
    AlbedoGlacierNirMin  = 0.4,           & ! Albedo of glacier in the NIR range at at the melting point
    AlbedoGlacierNirMax  = 0.75,          & ! Albedo of glacier in the NIR range at hard frost
    MinRadiation         = 1.0e-09,       & ! Minimum amount of radiation for albedo calculations [W/m2]
    ! Parameters for function vegetation_albedo_simple echam5
    AlbedoSnowCan        = 0.2              ! Albedo of snow on canopy
    !AlbedoSnowMin        = 0.4,   & ! Minimum possible snow albedo   ! Temporary, will depend on vegetation type
    !AlbedoSnowMax        = 0.8      ! Maximum possible snow albedo   ! Temporary

#endif
END MODULE mo_rad_constants
