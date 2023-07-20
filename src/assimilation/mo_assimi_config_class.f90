!> Contains structures and methods for assimi config
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
MODULE mo_assimi_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_math_constants,ONLY: eps4
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_assimi_config

  ! ======================================================================================================= !
  !> ASSIMI_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_assimi_config
    !> common parameters used by all models
    !>
    INTEGER                 :: ncanopy            !< number of canopy layers
    ! ----------------------------------------------------------------------------------------------------- !
    !> MODEL_JSBACH
    !>
    LOGICAL                 :: init_running_means ! initialize npp buffer for NLCC process
    REAL(wp), ALLOCATABLE   :: canopy_bound_lai(:)
    ! ----------------------------------------------------------------------------------------------------- !
    !> MODEL_QUINCY
    !>
    LOGICAL                 :: flag_optimal_Nfraction     !< optimise leaf internal N allocation
    LOGICAL                 :: flag_t_resp_acclimation    !< whether or not to use the respiration temperature acclimation factor
    LOGICAL                 :: flag_t_jmax_acclimation    !< acclimation of optimum temperature for Jmax
    CHARACTER(15)           :: canopy_layer_scheme        !< select scheme to calculate canopy layer thicknesses [standard|fapar]
    CHARACTER(15)           :: canopy_conductance_scheme  !< select scheme to calculate canopy conductance [medlyn|ballberry]
  CONTAINS
    PROCEDURE :: Init => Init_assimi_config
  END type t_assimi_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_config_class'

CONTAINS

  ! ======================================================================================================= !
  !> read ASSIMI_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_assimi_config(config)

  ! ======================================================================================================= !
  !> read ASSIMI_ namelist and init configuration parameters
  !>
  SUBROUTINE Init_assimi_config(config)
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_assimi_config), INTENT(inout) :: config
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    LOGICAL                     :: init_running_means
    INTEGER                     :: ncanopy
    INTEGER                     :: icanopy
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_canopy

    NAMELIST /jsb_assimi_nml/  &
      & active,                &
      & ic_filename,           &
      & bc_filename,           &
      & init_running_means

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_assimi_config'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Starting assimi configuration')
    ! ----------------------------------------------------------------------------------------------------- !

    ! Set defaults
    ! Set namelist defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_assimi.nc'
    ic_filename         = 'ic_land_assimi.nc'
    init_running_means  = .FALSE.
    ncanopy             = 3

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
#ifdef __QUINCY_STANDALONE__
    nml_unit    = position_nml('assimilation_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, assimilation_ctl)
#else
    nml_unit    = position_nml('jsb_assimi_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_assimi_nml)
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                       = active
    config%ic_filename                  = ic_filename
    config%bc_filename                  = bc_filename
    config%ncanopy                      = ncanopy_q_
    !! MODEL_JSBACH
    config%init_running_means           = init_running_means
    !! MODEL_QUINCY
    config%flag_optimal_Nfraction       = flag_optimal_Nfraction
    config%flag_t_resp_acclimation      = flag_t_resp_acclimation
    config%flag_t_jmax_acclimation      = flag_t_jmax_acclimation
    config%canopy_layer_scheme          = TRIM(canopy_layer_scheme)
    config%canopy_conductance_scheme    = TRIM(canopy_conductance_scheme)
    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%init_running_means  = init_running_means

    ! ncanopy is so far not a parameter in the assimi namelist but fixed to 3 above
    config%ncanopy = ncanopy

    ! Create vertical axis for canopy layers
    vgrid_canopy  => new_vgrid('canopy_layer', ZAXIS_GENERIC, config%ncanopy, &
      & longname = 'Canopy layers from ASSIMI_ JSBACH', &
      & levels = (/ (REAL(icanopy, kind=wp), icanopy=1, config%ncanopy) /), units = '')
    CALL register_vgrid(vgrid_canopy)

    ! Create vertical axis for canopy layers
    vgrid_canopy  => new_vgrid('canopy_layer', ZAXIS_GENERIC, config%ncanopy, &
      & longname = 'Canopy layers from ASSIMI_ JSBACH', &
      & levels = (/ (REAL(i, kind=wp), i=1, config%ncanopy) /), units = '')
    CALL register_vgrid(vgrid_canopy)
    
    ! Prepare canopy boundaries of lai_cl
    ALLOCATE(config%canopy_bound_lai(0:config%ncanopy))
    DO icanopy = 0,config%ncanopy
      config%canopy_bound_lai(icanopy) = REAL(icanopy, wp) / REAL(config%ncanopy, wp)
    END DO

    !! quincy
    !< canopy layer depth calculation
    !! 
    !! here the thickness of each canopy layer is calculated (saved in the vgrid as dz(:)) \n
    !! this thickness is not directly given (by e.g. using the unit meter), but represented by the leaf mass in this layer \n
    !! this leaf mass is measured as LAI [m2 m-2], that is, the actual thickness (in m) of canopy layer 2 depends on how much
    !!   vertical space (hight) needs to be covered to get the LAI specific to canopy layer 2
    !!
    !! layer 1 is the top canopy layer (reversed order of vgrid axis) \n
    !! and the value of the upper bound of a layer is larger than the value of the lower bound of the same layer
    !!
    ALLOCATE(canopy_layer_thickness_profile(config%ncanopy))
    SELECT CASE(TRIM(config%canopy_layer_scheme))
    CASE ("standard")
      canopy_layer_thickness_profile(:) = 0.5_wp
    CASE ("fapar")
      ! alternative canopy layer depth calculation, approximately proportional to FAPAR
      ! to reduce impact of actual canopy layering on PS calculation
      dAPAR = 1.0_wp / REAL(config%ncanopy, KIND = wp)
      DO icanopy = 1,config%ncanopy
        canopy_layer_thickness_profile(icanopy) = -1.0_wp / k_canopy_layer * LOG(1._wp - (dAPAR * icanopy) + eps4)
      END DO
    CASE DEFAULT
      WRITE(message_text,'(2a)') 'invalid canopy_layer_scheme', config%canopy_layer_scheme
      CALL finish(routine, message_text)
    END SELECT

    !> check validity of canopy-conductance scheme
    !>
    SELECT CASE(TRIM(config%canopy_conductance_scheme))
    CASE ("medlyn")
      CALL message(TRIM(routine), 'Using the canopy_conductance_scheme: '//TRIM(config%canopy_conductance_scheme)// &
        &                         ' the local lctlib_g1 variable is set to lctlib%g1_medlyn')
    CASE ("ballberry")
      CALL message(TRIM(routine), 'Using the canopy_conductance_scheme: '//TRIM(config%canopy_conductance_scheme)// &
        &                         ' the local lctlib_g1 variable is set to lctlib%g1_bberry')
    CASE DEFAULT
      WRITE(message_text,'(2a)') 'invalid canopy_conductance_scheme: ', TRIM(config%canopy_conductance_scheme)
      CALL finish(routine, message_text)
    END SELECT

    !> Create vertical axis - canopy layers
    !>
    !! these are "infrastructure" values, constant across PFTs and sites, and they do not interfere with tree height \n
    !! vgrid_canopy% lbounds, ubounds, and levels is calculated by 'new_vgrid()'
    !!
    !! note: the axis of the canopy_layer vgrid is used as "reversed", i.e., layer 1 is the top canopy layer \n
    !! levels is defined as the depth at the center of the layer: levels(:) = 0.5_wp * (lbounds(:) + ubounds(:))
    vgrid_canopy  => new_vgrid('canopy_layer', ZAXIS_GENERIC, config%ncanopy, &
      & longname = 'Canopy layers (leaf mass per layer as LAI) from ASSIMI_ QUINCY', &
      & units = 'm2 m-2 LAI', &
      ! levels=   , &
      ! lbounds=  , &
      ! ubounds=  , &
      & dz = canopy_layer_thickness_profile(:))
    CALL register_vgrid(vgrid_canopy)
  END SUBROUTINE Init_assimi_config_quincy

#endif
END MODULE mo_assimi_config_class
