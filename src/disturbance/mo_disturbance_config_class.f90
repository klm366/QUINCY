!>
!! @brief Contains structures and methods for disturbance config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (06.04.2017)
!! 2018-04-25, StW: Splitted configuration namelist in one general and a specific for each disturbance algorithm
!!
MODULE mo_disturb_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config
  USE mo_disturbance_constants, ONLY: windbreak_names, WINDBREAK_NONE, WINDBREAK_JSB3,               &
                                           fire_names,      FIRE_NONE,      FIRE_JSB3, FIRE_SPITFIRE

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_disturb_config

  TYPE, EXTENDS(t_jsb_config) :: t_disturb_config
    ! General settings
    INTEGER           :: fire_algorithm      ! Select fire algorithm: 0: none, 1: jsbach (default), 2: spitfire (Thonicke et al.(2001))
    INTEGER           :: windbreak_algorithm ! Select windbreak algorithm: 0: none, 1: jsbach (default)
    LOGICAL           :: lburn_pasture       ! Allow pasture pfts to burn
    ! Settings for windbreak-jsbach
    REAL(wp)          :: wnd_jsb3_threshold, &    ! factor by which the previous day maximum wind speed
                                                  ! must be larger than the climatological daily maximum
                                                  ! wind speed to allow any wind damage
                         wnd_jsb3_damage_scale    ! scaling factor for wind damage
    ! Settings for fire-jsbach
    REAL(wp)          :: fire_jsb3_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                         fire_jsb3_rel_hum_threshold, & ! maximal relative humidity for fire
                         fire_jsb3_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                         fire_jsb3_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                         fire_jsb3_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                         fire_jsb3_tau_grass            ! return period of fire for grass PFT [year] at 0% relative humidity
    ! Settings for fire-spitfire
    REAL(wp)          :: fire_spitfire_a_nd,              & ! Propensity of humans to ignite a fire
                         fire_spitfire_human_para,        & ! Parameter to tune the human ignitions
                         fire_spitfire_ign_para,          & ! Parameter to tune the sum of ignitions
                         fire_spitfire_SurfArea2Vol(3),   & ! Surface to volume ratio for 1hr 10hr and 100hr fuels respectively
                         fire_spitfire_moisture_scaling,  & ! Factor to convert SAtVR to fuel drying parameter
                         fire_spitfire_wind_max,          & ! Maximum effective windspeed
                         fire_spitfire_wind_slope,        & ! Determines the reduction of high windspeed
                         fire_spitfire_flaming2total,     & ! Factor to calculate total burning time from flaming phase duration
                         fire_spitfire_fract_green_active,& ! Fraction of leaves on woody types acting as fuel
                         fire_spitfire_const_ignitions,   & ! Set ignitions to this constant (-1 = no, default -1)
                         fire_spitfire_const_fract_fire,  & ! Set burned_frac to this constant (-1 = no, default -1)
                         fire_spitfire_fract_N_2_atm        ! Fraction of burned Nitrogen entering the atmos (rest SMINN)
    LOGICAL           :: fire_spitfire_llight_ground,     & ! File contains cloud to ground flashes only?
                         fire_spitfire_lwind_speed_limit, & ! Wind speed use for fire spreading is limited?
                         fire_spitfire_read_fuel_frac,    & ! Read fractions of fuel from file during initialization
                         fire_spitfire_lduration_popd,    & ! Fire duration depends on population density
                         fire_spitfire_lmortalityEQ1,     & ! Set mortality to 1 in all cases (default: .false.)
                         fire_spitfire_ldiag                ! Produce additional output (StW: General disturbance parameter?)
    ! Saved settings for fire-spitfire (not namelist controlled)
    CHARACTER(len=16) :: fire_spitfire_master_tile          ! At which tile is auxillary data stored?

    CONTAINS
      PROCEDURE :: Init => Init_disturb_config
  END type t_disturb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_config_class'
  

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize disturb process
  !!
  !! @param[inout]     config     Configuration type of process (t_disturb_config)
  !!
  SUBROUTINE Init_disturb_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_input_strings,      ONLY: GetStringIndex
    USE mo_input,              ONLY: InputGetLogicalFromString

    CLASS(t_disturb_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    CHARACTER(len=16)           :: fire_algorithm
    CHARACTER(len=16)           :: windbreak_algorithm
    CHARACTER(len=16)           :: lburn_pasture
 
    NAMELIST /jsb_disturb_nml/     &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         fire_algorithm,           &
         windbreak_algorithm,      &
         lburn_pasture

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_disturb_config'

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_disturb.nc'
    ic_filename            = 'ic_land_disturb.nc'
    fire_algorithm         = 'jsb3'
    windbreak_algorithm    = 'jsb3'
    lburn_pasture          = ''

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_disturb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) THEN
      READ(nml_unit, jsb_disturb_nml)
      CALL message(TRIM(routine), 'Obtained disturb configuration from: '//TRIM(config%namelist_filename))
    ENDIF

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%fire_algorithm         = GetStringIndex(TRIM(     fire_algorithm),     fire_names,ErrMsg=     'fire algorithm name') - 1
    config%windbreak_algorithm    = GetStringIndex(TRIM(windbreak_algorithm),windbreak_names,ErrMsg='windbreak algorithm name') - 1
    IF (LEN_TRIM(lburn_pasture) > 0) THEN
      config%lburn_pasture = InputGetLogicalFromString(TRIM(lburn_pasture),'lburn_pasture')
    ENDIF

    ! Read windbreak algorithm specific namelists here
    SELECT CASE (config%windbreak_algorithm)
      CASE (WINDBREAK_NONE)
        ! By definition nothing to do
      CASE (WINDBREAK_JSB3)
        CALL Init_disturb_windbreak_jsb3_config(config,nml_handler)
    END SELECT
    
    ! Read fire algorithm specific namelists here
    SELECT CASE (config%fire_algorithm)
      CASE (FIRE_NONE)
        ! By definition nothing to do
      CASE (FIRE_JSB3)
        CALL Init_disturb_fire_jsb3_config(config,nml_handler,lburn_pasture)
      CASE (FIRE_SPITFIRE)
        CALL Init_disturb_fire_spitfire_config(config,nml_handler,lburn_pasture)
    END SELECT
    
    CALL close_nml(nml_handler)

  END SUBROUTINE Init_disturb_config

  ! Initialises configuration of windbreak algorithm jsb3
  SUBROUTINE Init_disturb_windbreak_jsb3_config(config,nml_handler)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_disturb_config), INTENT(inout) :: config
    INTEGER,                 INTENT(in)    :: nml_handler

    ! windbreak-jsb3 namelist parameters:
    REAL(wp) :: threshold, & ! factor by which the previous day maximum wind speed
                             ! must be larger than the climatological daily maximum
                             ! wind speed to allow any wind damage
                damage_scale ! scaling factor for wind damage

    NAMELIST /jsb_disturb_windbreak_jsb3_nml/ &
         threshold,                           &
         damage_scale

    INTEGER :: nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_disturb_windbreak_jsb3_config'

    ! Set default values
    threshold    =  2.25_wp   ! factor by which the previous day maximum wind speed
                              ! must be larger than the climatological daily maximum
                              ! wind speed to allow any wind damage
    damage_scale =  5.e-03_wp ! scaling factor for wind damage

    ! Read namelist
    nml_unit = position_nml('jsb_disturb_windbreak_jsb3_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) THEN
      READ(nml_unit, jsb_disturb_windbreak_jsb3_nml)
      CALL message(TRIM(routine), 'Obtained windbreak-jsb3 configuration from: '//TRIM(config%namelist_filename))
    ENDIF

    config%wnd_jsb3_threshold    = threshold
    config%wnd_jsb3_damage_scale = damage_scale

  END SUBROUTINE Init_disturb_windbreak_jsb3_config

  ! Initialises configuration of fire algorithm jsb3
  SUBROUTINE Init_disturb_fire_jsb3_config(config,nml_handler,lburn_pasture)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_disturb_config), INTENT(inout) :: config
    INTEGER,                 INTENT(in)    :: nml_handler
    CHARACTER(len=*),        INTENT(in)    :: lburn_pasture

    ! fire-jsb3 namelist parameters:
    REAL(wp)                    :: litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                                   rel_hum_threshold, & ! maximal relative humidity for fire
                                   minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                                   minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                                   tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                                   tau_grass            ! return period of fire for grass PFT [year] at 0% relative humidity

    NAMELIST /jsb_disturb_fire_jsb3_nml/ &
         litter_threshold,               &
         rel_hum_threshold,              &
         minimum_woody,                  &
         minimum_grass,                  &
         tau_woody,                      &
         tau_grass

    INTEGER :: nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_disturb_fire_jsb3_config'

    ! Set default values
    litter_threshold  = 16.67_wp   ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
    rel_hum_threshold = 70._wp     ! maximal relative humidity for fire
    minimum_woody     =  0.002_wp  ! minimal fraction of act_fpc of woody PFT to be burned each year
    minimum_grass     =  0.006_wp  ! minimal fraction of act_fpc of grass PFT to be burned each year
    tau_woody         =  6.0_wp    ! return period of fire for woody PFT [year] at 0% relative humidity
    tau_grass         =  2.0_wp    ! return period of fire for grass PFT [year] at 0% relative humidity
    IF (LEN_TRIM(lburn_pasture)==0) THEN
      config%lburn_pasture = .FALSE.
    ENDIF

    ! Read namelist
    nml_unit = position_nml('jsb_disturb_fire_jsb3_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) THEN
      READ(nml_unit, jsb_disturb_fire_jsb3_nml)
      CALL message(TRIM(routine), 'Obtained fire-jsb3 configuration from: '//TRIM(config%namelist_filename))
    ENDIF

    config%fire_jsb3_litter_threshold  = litter_threshold
    config%fire_jsb3_rel_hum_threshold = rel_hum_threshold
    config%fire_jsb3_minimum_woody     = minimum_woody
    config%fire_jsb3_minimum_grass     = minimum_grass
    config%fire_jsb3_tau_woody         = tau_woody
    config%fire_jsb3_tau_grass         = tau_grass

  END SUBROUTINE Init_disturb_fire_jsb3_config

  ! Initialises configuration of fire algorithm spitfire
  SUBROUTINE Init_disturb_fire_spitfire_config(config,nml_handler,lburn_pasture)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_disturb_config), INTENT(inout) :: config
    INTEGER,                 INTENT(in)    :: nml_handler
    CHARACTER(len=*),        INTENT(in)    :: lburn_pasture

    ! Local namelist parameters
    REAL(wp) :: a_nd
    REAL(wp) :: human_para
    REAL(wp) :: ign_para
    REAL(wp) :: wind_max
    REAL(wp) :: wind_slope
    REAL(wp) :: flaming2total
    REAL(wp) :: SurfArea2Vol(3)
    REAL(wp) :: moisture_scaling
    REAL(wp) :: fract_green_active
    REAL(wp) :: const_ignitions
    REAL(wp) :: const_fract_fire
    REAL(wp) :: fract_N_2_atm
    LOGICAL  :: llight_ground
    LOGICAL  :: read_fuel_frac
    LOGICAL  :: lduration_popd
    LOGICAL  :: lmortalityEQ1
    LOGICAL  :: ldiag
    LOGICAL  :: lwind_speed_limit

    NAMELIST /jsb_disturb_fire_spitfire_nml/ &
      llight_ground,                         &
      lwind_speed_limit,                     &
      a_nd,                                  &
      human_para,                            &
      ign_para,                              &
      SurfArea2Vol,                          &
      moisture_scaling,                      &
      wind_max,                              &
      wind_slope,                            &
      read_fuel_frac,                        &
      lduration_popd,                        &
      flaming2total,                         &
      lmortalityEQ1,                         &
      const_ignitions,                       &
      const_fract_fire,                      &
      fract_green_active,                    &
      fract_N_2_atm,                         &
      ldiag

    INTEGER :: nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_disturb_fire_jsb3_config'

    ! Set default values
    SurfArea2Vol(:)   = (/ 66.0_wp,3.58_wp,0.98_wp /)
    a_nd              =     0.251_wp
    human_para        =     0.4_wp
    ign_para          =     0.22_wp
    wind_max          =   150._wp
    wind_slope        =     1.5_wp
    moisture_scaling  = 30000.0_wp
    flaming2total     = 5.0_wp
    fract_N_2_atm     = 0.5_wp
    llight_ground     = .FALSE.
    lwind_speed_limit = .TRUE.
    lduration_popd    = .TRUE.
    fract_green_active= 0.3_wp
    lmortalityEQ1     = .FALSE.
    const_ignitions   = -1.0_wp
    const_fract_fire  = -1.0_wp
    ldiag             = .FALSE.
    IF (LEN_TRIM(lburn_pasture)==0) THEN
      config%lburn_pasture = .TRUE.
    ENDIF
    config%fire_spitfire_master_tile = ''

    ! Read namelist
    nml_unit = position_nml('jsb_disturb_fire_spitfire_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) THEN
      READ(nml_unit, jsb_disturb_fire_spitfire_nml)
      CALL message(TRIM(routine), 'Obtained fire-spitfire configuration from: '//TRIM(config%namelist_filename))
    ENDIF

    config%fire_spitfire_a_nd              = a_nd
    config%fire_spitfire_human_para        = human_para
    config%fire_spitfire_ign_para          = ign_para
    config%fire_spitfire_wind_max          = wind_max
    config%fire_spitfire_wind_slope        = wind_slope
    config%fire_spitfire_flaming2total     = flaming2total
    config%fire_spitfire_SurfArea2Vol(:)   = SurfArea2Vol(:)
    config%fire_spitfire_moisture_scaling  = moisture_scaling
    config%fire_spitfire_fract_green_active= fract_green_active
    config%fire_spitfire_const_ignitions   = const_ignitions
    config%fire_spitfire_const_fract_fire  = const_fract_fire
    config%fire_spitfire_fract_N_2_atm     = fract_N_2_atm
    config%fire_spitfire_llight_ground     = llight_ground
    config%fire_spitfire_read_fuel_frac    = read_fuel_frac
    config%fire_spitfire_lwind_speed_limit = lwind_speed_limit
    config%fire_spitfire_lduration_popd    = lduration_popd
    config%fire_spitfire_lmortalityEQ1     = lmortalityEQ1
    config%fire_spitfire_ldiag             = ldiag

  END SUBROUTINE Init_disturb_fire_spitfire_config

#endif
END MODULE mo_disturb_config_class
