!>
!! @brief Contains structures and methods for config of the assimilation process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! integration into SLM by Jan Engel, MPI BGC (20160722)
!!
MODULE mo_assimi_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp, eps4
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_assimi_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the assimilation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_assimi_config
    LOGICAL                 :: flag_optimal_Nfraction     !< optimise leaf internal N allocation
    LOGICAL                 :: flag_t_resp_acclimation    !< whether or not to use the respiration temperature acclimation factor
    LOGICAL                 :: flag_t_jmax_acclimation    !< acclimation of optimum temperature for Jmax 
    INTEGER                 :: ncanopy                    !< number of canopy layer simulated by the model
    CHARACTER(15)           :: canopy_layer_scheme        !< select scheme to calculate depth of each single canopy layer
                                                          ! "standard|fapar"           
    REAL(wp), ALLOCATABLE   :: canopy_layer_depth(:)      !< depth of each single canopy layer
  CONTAINS
    PROCEDURE :: Init => Init_assimi_config
  END type t_assimi_config


  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_assimi_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_assimi_config(config)   ! 'config' == t_assimi_config, see ABSTRACT INTERFACE -> SUBROUTINE Init_config(config)

!    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
!    using the old JSBACH3 functions instead of wrappers implemented in JSBACH4 -> mo_jsb_namelist_iface
    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    ! 0.1 InOut
    CLASS(t_assimi_config), INTENT(inout) :: config    !< config type for assimi

    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: flag_optimal_Nfraction     , &
                                   flag_t_resp_acclimation    , &
                                   flag_t_jmax_acclimation
    CHARACTER(15)               :: canopy_layer_scheme
    INTEGER                     :: ncanopy    , &
                                   icanopy
    REAL(wp)                    :: dAPAR
    REAL(wp)                    :: k_canopy_layer                !< extinction coefficient to calculate canopy layer depth
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_assimi_config'

    NAMELIST /assimilation_ctl/    &
        active,                    &
        ic_filename,               &
        bc_filename,               &
        flag_optimal_Nfraction,    &
        flag_t_resp_acclimation,   &
        flag_t_jmax_acclimation,   &
        ncanopy,                   &
        canopy_layer_scheme

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    ! extinction coefficient for canopy layer depth calculation
    k_canopy_layer = 0.7_wp

    IF (debug_on()) CALL message(TRIM(routine), 'Starting assimi configuration')

    ! Set defaults
    active                         = .TRUE.
    bc_filename                    = 'bc_land_assimi.nc'
    ic_filename                    = 'ic_land_assimi.nc'
    flag_optimal_Nfraction         = .FALSE.      ! .TRUE. .FALSE.
    flag_t_resp_acclimation        = .FALSE.      ! .TRUE. .FALSE.
    flag_t_jmax_acclimation        = .FALSE.      ! .TRUE. .FALSE.
    ncanopy                        = 20           ! (default for FAPAR based canopy layers would be 10)
    canopy_layer_scheme            = "standard"   ! 

    ! read the namelist
    ! &ASSIMILATION_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('ASSIMILATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, assimilation_ctl)

    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%active                        = active
    config%ic_filename                   = ic_filename
    config%bc_filename                   = bc_filename
    config%flag_optimal_Nfraction        = flag_optimal_Nfraction
    config%flag_t_resp_acclimation       = flag_t_resp_acclimation
    config%flag_t_jmax_acclimation       = flag_t_jmax_acclimation
    config%ncanopy                       = ncanopy
    config%canopy_layer_scheme           = TRIM(canopy_layer_scheme)


    !< canopy layer depth calculation
    !! 
    ALLOCATE(config%canopy_layer_depth(ncanopy))

    SELECT CASE(TRIM(config%canopy_layer_scheme))
      CASE ("standard")
        config%canopy_layer_depth(:) = 0.5_wp

      CASE ("fapar")
        ! alternative canopy layer depth calculation, approximately proportional to FAPAR to  
        ! reduce impact of actual canopy layering on PS calculation
        dAPAR = 1._wp/REAL(ncanopy,KIND=wp)
        DO icanopy=1,ncanopy
          config%canopy_layer_depth(icanopy) = -1.0_wp/k_canopy_layer * LOG(1._wp-(dAPAR*icanopy)+eps4)
        END DO

      CASE DEFAULT
        WRITE(message_text,'(2a)') 'invalid canopy_layer_scheme',canopy_layer_scheme
        CALL finish("mo_assimi_config_class:Init_assimi_config",message_text) 
      
    END SELECT

  END SUBROUTINE Init_assimi_config

#endif
END MODULE mo_assimi_config_class
