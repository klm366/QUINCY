!>
!! @brief Contains structures and methods for disturbance config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (06.04.2017)
!!
MODULE mo_disturb_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

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
    REAL(wp)          :: wnd_threshold, &    ! factor by which the previous day maximum wind speed
                                             ! must be larger than the climatological daily maximum
                                             ! wind speed to allow any wind damage
                         wnd_damage_scale    ! scaling factor for wind damage
    ! Settings for fire-jsbach
    REAL(wp)          :: fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                         fire_rel_hum_threshold, & ! maximal relative humidity for fire
                         fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                         fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                         fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                         fire_tau_grass            ! return period of fire for grass PFT [year] at 0% relative humidity

    CONTAINS
      PROCEDURE :: Init => Init_disturb_config
  END type t_disturb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_config_class'
  
  CHARACTER(len=*), PARAMETER :: windbreak_names = 'none jsb3'
  INTEGER, PUBLIC,  PARAMETER :: WINDBREAK_NONE  = 0
  INTEGER, PUBLIC,  PARAMETER :: WINDBREAK_JSB3  = 1

  CHARACTER(len=*), PARAMETER :: fire_names      = 'none jsb3 spitfire'
  INTEGER, PUBLIC,  PARAMETER :: FIRE_NONE       = 0
  INTEGER, PUBLIC,  PARAMETER :: FIRE_JSB3       = 1
  INTEGER, PUBLIC,  PARAMETER :: FIRE_SPITFIRE   = 2


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
!        CALL Init_disturb_fire_spitfire_config(config,nml_handler,lburn_pasture)
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

    config%wnd_threshold    = threshold
    config%wnd_damage_scale = damage_scale

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

    config%fire_litter_threshold  = litter_threshold
    config%fire_rel_hum_threshold = rel_hum_threshold
    config%fire_minimum_woody     = minimum_woody
    config%fire_minimum_grass     = minimum_grass
    config%fire_tau_woody         = tau_woody
    config%fire_tau_grass         = tau_grass

  END SUBROUTINE Init_disturb_fire_jsb3_config

#endif
END MODULE mo_disturb_config_class
