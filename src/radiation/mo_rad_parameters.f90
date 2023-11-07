!>
!! @brief Contains parameters for the radiation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-15)
!! copy all parameters to here from mo_rad_constantss by Jan Engel (2023-04-05)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_rad_parameters
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters

  IMPLICIT NONE
  PUBLIC

  ! ======================================================================================================= !
  !> 1.0 general RAD_ parameters
  !>
  REAL(wp), SAVE :: &
   & MinRadiation     = def_parameters      !< Minimum amount of radiation for albedo calculations [W/m2]

  ! ======================================================================================================= !
  !> 2.0 JSBACH specific parameters
  !>
  REAL(wp), SAVE :: &
   & Epar                   = def_parameters,    & ! Energy content of PAR [J / mol(photons)]=(4.6 mol/MJ PAR)**-1.  R: for update_par
   & SoilReflectivityParMin = def_parameters,      & ! Minimum soil reflectivity in PAR region.  R: for update_par
  !  & LaiMax = def_parameters,                       & ! Maximum LAI (used for nitrogen scaling)
   & LaiMin = def_parameters,                    & ! Min Lai in PAR computation for calc_fapar_lai_cl
   & LaiLimit = def_parameters,                     & ! Min LAI: used for nitro scaling and estimation of fract cover for calc_fapar_lai_cl
   & FcMax = def_parameters,                       & ! Max fractional vegetation cover. R: for calc_fapar_lai_cl
   & FcMin = def_parameters,                     & ! Min fractional vegetation cover. R: for calc_fapar_lai_cl
   & ZenithMinPar = def_parameters,              & ! Min cos of zenith angle for which PAR is calculated. R: for calc_fapar_lai_cl
   !     ZenithMin = def_parameters,            & ! Check for solar zenith angle > 89 degrees
   !     SoilReflectivityParMin = def_parameters,     & ! Minimum soil reflectivity in PAR region
   !
   !     & minStomaConductance = def_parameters,        & ! Minimum stomatal conductance [mol H2O /(m^2 s) ??]
   ! Parameters for spatially constant albedo of the soil surface (UseAlbedoSoilConst)
   & AlbedoSoilConstVis = def_parameters,          & ! Albedo constant alternatively usable for alb_vis_mineralsoil
   & AlbedoSoilConstNir = def_parameters,          & ! See above, for the nir range
                                !
   & SkyViewFactor      = def_parameters,          & ! Constant in calculating the sky view fractions (only to scale the e-function)
   & ZenithAngleFactor  = def_parameters,           & ! Factor in solar zenith angle dependence of snow albedo
                                !                                 ! (the increase of snow albedo is the higher this factor
                                ! Parameters for albedo
   & AlbedoInitial    = def_parameters,            & ! initial value for albedo in the whole solar range
   & AlbedoVisInitial = def_parameters,            & ! initial value for albedo in the visible range
   & AlbedoNirInitial = def_parameters,            & ! initial value for albedo in the near infrared
                                ! Parameters of snow age scheme albedo
   & AlbedoSnowVisMax_age   = def_parameters,     & ! Maximum albedo of fresh snow in the visible range
   & AlbedoSnowNirMax_age   = def_parameters,     & ! Maximum albedo of fresh snow in the NIR range
   & AlbedoSnowVisAge       = def_parameters,     & ! Maximal rel. reduction of snow albedo by aging in the visible range
   & AlbedoSnowNirAge       = def_parameters,      & ! Maximal rel. reduction of snow albedo by aging in the NIR range
   & AlbedoCanopySnow_age   = def_parameters,     & ! Albedo of snow covered canopy for snow age scheme (BATS)
   & AlbedoSnowAngle        = def_parameters,      & ! Maximal rel. reduction of snow absorption by large solar zenith angle
                                !
                                ! Parameters of snow temperature scheme albedo
   & AlbedoSnowVisMax_temp   = def_parameters,    & ! Maximum albedo of fresh snow in the visible range
   & AlbedoSnowNirMax_temp   = def_parameters,    & ! Maximum albedo of fresh snow in the NIR range
   & AlbedoSnowVisMin_temp   = def_parameters,    & ! Minimum albedo of fresh snow in the visible range
   & AlbedoSnowNirMin_temp   = def_parameters,     & ! Minimum albedo of fresh snow in the NIR range
   & TempAlbedoSnowMax       = def_parameters,     & ! Maximum snow albedo at this temperature below melting point of H2O
   & AlbedoCanopySnow_temp   = def_parameters,    & ! Albedo of snow covered canopy for temperature scheme
    !
    ! Parameters for lake albedo
   & AlbedoLakeWater         = def_parameters,    & ! Albedo of open lake water
    ! @todo: These should be set dependent on resolution (see mo_surface_ice:init_albedo_ice in ECHAM6.2)
    ! The values below correspond to the T63 settings
   & AlbedoLakeIceMin        = def_parameters,    &
   & AlbedoLakeIceMax        = def_parameters,    &
   & AlbedoLakeSnowMin       = def_parameters,    &
   & AlbedoLakeSnowMax       = def_parameters,    &
    !
    ! Parameters of glacier scheme albedo
   & AlbedoGlacierVisMin  = def_parameters,       & ! Albedo of glacier in the visible range at the melting point
   & AlbedoGlacierVisMax  = def_parameters,       & ! Albedo of glacier in the visible range at hard frost
   & TempAlbedoGlacierMax = def_parameters,        & ! Maximum glacier albedo at this temperature below melting point of H2O
   & AlbedoGlacierNirMin  = def_parameters,       & ! Albedo of glacier in the NIR range at at the melting point
   & AlbedoGlacierNirMax  = def_parameters,       & ! Albedo of glacier in the NIR range at hard frost
    ! Parameters for function vegetation_albedo_simple echam5
   & AlbedoCanopySnow_simple = def_parameters,     & ! Albedo of snow on canopy
   & AlbedoSnowMin_simple    = def_parameters,     & ! Minimum possible snow albedo
   & AlbedoSnowMax_simple    = def_parameters        ! Maximum possible snow albedo

  ! ======================================================================================================= !
  !> 3.0 QUINCY specific parameters
  !>
  ! Parameters required for albedo and fapar calculation
  REAL(wp), SAVE :: & 
    & kbl0_vis = def_parameters, &               !< extinction coefficient over black leaves 
    & kdf0_vis = def_parameters, &               !< extinction coefficient for diffuse PAR (Spitters eq. 7) 
    & kbl0_nir = def_parameters, &               !< extinction coefficient over black leaves in NIR range 
    & kdf0_nir = def_parameters, &               !< extinction coefficient for diffuse NIR (Spitters eq. 7) 
    & rho2sbeta = def_parameters, &              !< scaling factor of solar angle in reflection calculation 
    & k_csf = def_parameters, &                  !< scaling factor in the crown-shape effect on clumping (Campbell & Norman, 1998, eq 15.35)
    & min_cos_zenith_angle_rad = def_parameters  !< minimum coszine zenith angle for radiation calculation

  REAL(wp), SAVE :: & 
    & def_alb_vis_soil = def_parameters, &  !< default VIS albedo for soil
    & def_alb_nir_soil = def_parameters, &  !< default NIR albedo for soil
    & def_alb_vis_snow = def_parameters, &  !< default VIS albedo for snow
    & def_alb_nir_snow = def_parameters     !< default NIR albedo for snow

  ! conversion factors
  REAL(wp), SAVE :: &
    & rad2ppfd = def_parameters             !< micro-mol /J

  ! red to far red ratio calculations
  REAL(wp), SAVE :: &
    & rfr_ratio_toc = def_parameters, &     !< top of the canopy red:farred ratio  (Kull & Kruij 1998) 
    & k_r2fr_chl = def_parameters           !< extinction of red:farred ratio on leaf chlorophyll (1/mmol N); 
                                            !! (Kull & Krujit 1998 for the far red value, tuned to get a R:FR ratio
                                            !! of 0.15 at full canopy closure (a common value) 

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_rad_parameters'

CONTAINS

  ! ======================================================================================================= !
  !> initialize QUINCY specific parameters for the process: radiation
  !>
  SUBROUTINE init_rad_parameters_jsbach
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_rad_parameters_jsbach'

    ! ...
    Epar                    = 2.2E5_wp
    SoilReflectivityParMin  = 0.0_wp
    LaiMin                  = 1.E-9_wp
    LaiLimit                = 3._wp
    FcMax                   = 1.0_wp
    FcMin                   = 1.E-3_wp
    ZenithMinPar            = 1.E-3_wp
    ! Parameters for spatially constant albedo of the soil surface (UseAlbedoSoilConst)
    AlbedoSoilConstVis      = 0.2_wp
    AlbedoSoilConstNir      = 0.2_wp
    SkyViewFactor           = 0.5_wp
    ZenithAngleFactor       = 2._wp
    ! Parameters for albedo
    AlbedoInitial           = 0.2_wp
    AlbedoVisInitial        = 0.1_wp
    AlbedoNirInitial        = 0.3_wp
    MinRadiation            = 1.0e-09_wp
    ! Parameters of snow age scheme albedo
    AlbedoSnowVisMax_age    = 0.90_wp
    AlbedoSnowNirMax_age    = 0.60_wp
    AlbedoSnowVisAge        = 0.15_wp
    AlbedoSnowNirAge        = 0.5_wp
    AlbedoCanopySnow_age    = 0.20_wp
    AlbedoSnowAngle         = 0.4_wp
    ! Parameters of snow temperature scheme albedo
    AlbedoSnowVisMax_temp   = 0.90_wp
    AlbedoSnowNirMax_temp   = 0.65_wp
    AlbedoSnowVisMin_temp   = 0.52_wp
    AlbedoSnowNirMin_temp   = 0.3_wp
    TempAlbedoSnowMax       = 5.0_wp
    AlbedoCanopySnow_temp   = 0.20_wp
    ! Parameters for lake albedo
    AlbedoLakeWater         = 0.07_wp
    ! The values below correspond to the T63 settings  @TODO - update to ICON grid ?
    AlbedoLakeIceMin        = 0.60_wp
    AlbedoLakeIceMax        = 0.75_wp
    AlbedoLakeSnowMin       = 0.65_wp
    AlbedoLakeSnowMax       = 0.80_wp
    ! Parameters of glacier scheme albedo
    AlbedoGlacierVisMin     = 0.76_wp
    AlbedoGlacierVisMax     = 0.89_wp
    TempAlbedoGlacierMax    = 5.0_wp
    AlbedoGlacierNirMin     = 0.34_wp
    AlbedoGlacierNirMax     = 0.72_wp
    ! Parameters for function vegetation_albedo_simple echam5
    AlbedoCanopySnow_simple = 0.2_wp
    AlbedoSnowMin_simple    = 0.4_wp
    AlbedoSnowMax_simple    = 0.8_wp
  END SUBROUTINE init_rad_parameters_jsbach

  ! ======================================================================================================= !
  !> initialize QUINCY specific parameters for the process: radiation
  !>
  SUBROUTINE init_rad_parameters_quincy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_rad_parameters_quincy'

    ! Parameters required for albedo and fapar calculation
    kbl0_vis                  = 0.5_wp        !< Spitters 1986 
    kdf0_vis                  = 0.8_wp        !< Spitters 1986 
    kbl0_nir                  = 0.5_wp        !< Spitters 1986 
    kdf0_nir                  = 0.8_wp        !< Spitters 1986 
    rho2sbeta                 = 1.6_wp        !< Spitters 1986 
    k_csf                     = 2.2_wp        !< Campbell & Norman 1998, eq 15.35
    min_cos_zenith_angle_rad  = 0.1736482_wp  !< corresponding to 10 deg solar elevation, below which the 
                                              !! radiation equations are no longer valid
    def_alb_vis_soil          = 0.15_wp       !< default from Bonan 2008
    def_alb_nir_soil          = 0.30_wp       !< default from Bonan 2008
    def_alb_vis_snow          = 0.6_wp        !< Literature research: from 0.2(old snow) to 0.9(new snow) | @TODO: could vary with snow age 
    def_alb_nir_snow          = 0.6_wp        !< Literature research: from 0.2(old snow) to 0.9(new snow) | @TODO: could vary with snow age
    ! conversion factors
    rad2ppfd                  = 4.6_wp        !< Monteith & Unsworth 1995 
    ! red to far red ratio calculations
    rfr_ratio_toc             = 1.2_wp                !< Kull & Kruijt 1998
    k_r2fr_chl                = 0.0012_wp - 0.012_wp  !< Kull & Kruijt 1998
    ! parameters also used in jsbach
    MinRadiation              = 1.0e-09_wp
  END SUBROUTINE init_rad_parameters_quincy

#endif
END MODULE mo_rad_parameters
