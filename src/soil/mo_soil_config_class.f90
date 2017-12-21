!>
!! @brief Contains structures and methods for soil config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! 2nd version by Jan Engel, MPI BGC (20160913) -> integration into SLM
!!
MODULE mo_soil_config_class
!#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_kind,              ONLY: eps8
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_soil_config
  PUBLIC :: max_soil_layers

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the soil process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_soil_config
    INTEGER               :: nsoil_energy         !< number of soil layers for soil energy calculations
    INTEGER               :: nsoil_water          !< number of soil layers for soil moisture calculations
    REAL(wp)              :: soil_depth           !< actual soil and rooting depth
    REAL(wp), ALLOCATABLE :: soil_depth_energy(:) !< soil depth of the layers (lower bound)
    REAL(wp), ALLOCATABLE :: soil_depth_water(:)
    REAL(wp), ALLOCATABLE :: dz_energy(:)         !< layer thicknesses
    REAL(wp), ALLOCATABLE :: dz_water(:)
    REAL(wp)              :: wsr_capacity         !< Water holding capacity of ground skin reservoir [m water equivalent]
    REAL(wp)              :: wsn_capacity         !< Water holding capacity of ground snow reservoir (max snow depth) [m water equivalent]
    REAL(wp)              :: soil_awc_prescribe, &
                             soil_theta_prescribe
    REAL(wp)              :: soil_sand            !< soil sand proportion - site specific from forcing info file
    REAL(wp)              :: soil_silt            !< soil silt proportion - site specific from forcing info file
    REAL(wp)              :: soil_clay            !< soil clay proportion - recalculated from: clay = 1.0 -sand -silt
    REAL(wp)              :: bulk_density    
    CHARACTER(15)         :: soil_moisture_scheme !< soil moisture "prescribed|bucke|layer"
  CONTAINS
    PROCEDURE :: Init => Init_soil_config
  END TYPE t_soil_config

  INTEGER, PARAMETER :: max_soil_layers = 20  ! consistent with JSBACH4

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_soil_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_soil_config(config)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    ! 0.1 InOut
    CLASS(t_soil_config), INTENT(inout) :: config    !< config type for soil

    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL  :: active
    INTEGER  :: nsoil_energy    , &
                nsoil_water
    REAL(wp) :: dz_energy(max_soil_layers)
    REAL(wp) :: dz_water(max_soil_layers)
    REAL(wp) :: soil_depth_profile(max_soil_layers)
    REAL(wp) :: wsr_capacity    , &
                wsn_capacity    , &
                soil_depth   
    REAL(wp) :: soil_awc_prescribe    , &
                soil_theta_prescribe
    REAL(wp) :: soil_sand       , &
                soil_silt       , &
                soil_clay
    REAL(wp) :: bulk_density
    CHARACTER(15) :: soil_moisture_scheme
    INTEGER  :: i   ! looping
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    REAL(wp) :: k_soil_profile  = 0.25_wp
    REAL(wp) :: min_layer_depth = 0.065_wp

    NAMELIST /soil_ctl/       & 
                          active          , &
                          ic_filename     , &
                          bc_filename     , &
                          nsoil_energy    , &
                          nsoil_water     , &
                          dz_energy       , &
                          dz_water        , &
                          wsr_capacity    , &
                          wsn_capacity    , &
                          soil_depth      , &
                          soil_awc_prescribe        , &
                          soil_theta_prescribe      , &
                          soil_sand       , &
                          soil_silt       , &
                          soil_clay       , &
                          bulk_density    , &
                          soil_moisture_scheme

    ! variables for reading model-options from namelist
    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_soil_config'

    CALL message(TRIM(routine), 'Starting soil configuration')

    ! Set defaults
    active                        = .TRUE.
    ic_filename                   = 'ic_land_soil.nc'
    bc_filename                   = 'bc_land_soil.nc'
    nsoil_energy                  = 5
    nsoil_water                   = 5
    dz_energy(1:5)                = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)
    dz_energy(6:max_soil_layers)  = 0._wp
    dz_water (1:5)                = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)
    dz_water (6:max_soil_layers)  = 0._wp 
    wsr_capacity                  = 2.E-4_wp
    wsn_capacity                  = 2.E-4_wp
    soil_depth                    = 5.7_wp 
    soil_awc_prescribe            = 300.0_wp
    soil_theta_prescribe          = 1.0_wp 
    soil_sand                     = 0.3_wp
    soil_silt                     = 0.4_wp
    soil_clay                     = 0.3_wp
    bulk_density                  = 1500.0_wp
    soil_moisture_scheme          = "bucket" ! prescribed bucket layer 

    ! read the namelist
    ! &SOIL_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('SOIL_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, soil_ctl)

    CALL close_nml(nml_handler)

    ! allocate vectors
    ALLOCATE(config%soil_depth_energy(nsoil_energy))
    ALLOCATE(config%soil_depth_water (nsoil_water ))
    ALLOCATE(config%dz_energy(nsoil_energy))
    ALLOCATE(config%dz_water (nsoil_water ))

    ! check initialisation of number of soil layers
    IF (nsoil_water .LT. 5 .OR. nsoil_water .GT. 20) & 
        CALL finish('Init_soil_config','Number of soil layers needs to be between 5 and 20.')
    IF (nsoil_water .NE. nsoil_energy) & 
        CALL finish('Init_soil_config','Initialization with different numbers of soil layers for energy and water not supported.')

    ! If not using 5 layers, estimate depth profile to mimick standard JSB4 profile
    IF(nsoil_water > 5) THEN
       DO i=1,nsoil_water
         soil_depth_profile(i) = min_layer_depth * & 
                                 exp(k_soil_profile * i * REAL(max_soil_layers)/REAL(nsoil_water)) - & 
                                 min_layer_depth
       ENDDO 
       dz_energy(1) = min_layer_depth
       DO i=2,nsoil_water
         dz_energy(i) = MAX(min_layer_depth,soil_depth_profile(i) - soil_depth_profile(i-1))
       ENDDO
       dz_water(:) = dz_energy(:)
    ENDIF

    IF (soil_depth .LE. 0.0_wp) CALL finish('Init_soil_config','Initialization with non positive soil depth')
    IF (soil_depth .GT. SUM(dz_energy)) CALL finish('Init_soil_config','Initialization of soil depth larger than soil profile not possible')

    ! check input for soil_moisture_scheme 
    SELECT CASE (TRIM(soil_moisture_scheme))
     CASE ("prescribed","bucket","layer")
     !CASE ("prescribed","bucket")
     !CASE ("layer") 
     !  WRITE(message_text,'(2a)') soil_moisture_scheme,'is not implemented yet, sorry....'
     !  CALL finish("mo_soil_config_class:init_soil_config",message_text)
     CASE DEFAULT
       WRITE(message_text,'(2a)') 'invalid soil_moisture_scheme',soil_moisture_scheme
       CALL finish("mo_soil_config_class:init_soil_config",message_text)
    ENDSELECT

    ! pass values as read from file
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%nsoil_energy              = nsoil_energy
    config%nsoil_water               = nsoil_water
    config%soil_awc_prescribe        = soil_awc_prescribe 
    config%soil_theta_prescribe      = soil_theta_prescribe 
    config%soil_moisture_scheme      = soil_moisture_scheme

    config%soil_depth_energy(1)      = dz_energy(1)
    DO i=2,nsoil_energy
      config%soil_depth_energy(i) = config%soil_depth_energy(i-1) + dz_energy(i)
    END DO
    config%soil_depth_water(1)       = dz_water(1)
    DO i=2,nsoil_water
      config%soil_depth_water(i) = config%soil_depth_water(i-1) + dz_water(i)
    END DO

    config%dz_energy(1:nsoil_energy) = dz_energy(1:nsoil_energy)
    config%dz_water (1:nsoil_water)  = dz_water(1:nsoil_water)
    config%wsr_capacity              = wsr_capacity
    config%wsn_capacity              = wsn_capacity
    config%soil_depth                = soil_depth 

    config%soil_sand                 = soil_sand
    config%soil_silt                 = soil_silt
    !config%soil_clay                 = soil_clay must be calculated from soil_sand & soil_silt to guarantee: 1=sand+silt+clay
    config%soil_clay                 = 1.0_wp - soil_sand - soil_silt
    IF(config%soil_clay < eps8) THEN
       WRITE(message_text,'(2a)') 'invalid proportion of soil_clay (value < eps8)',routine
       CALL finish("mo_soil_config_class:init_soil_config",message_text)
    END IF
    config%bulk_density              = bulk_density

  END SUBROUTINE Init_soil_config

!#endif
END MODULE mo_soil_config_class
