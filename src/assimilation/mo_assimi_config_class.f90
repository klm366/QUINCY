!> Summary: JSBACH & QUINCY assimilation process config
!>
!>#### ICON-Land
!>
!> ------------------------------------------------------------------------
!> Copyright (C) 2004-2024, MPI-M, MPI-BGC, DWD, DKRZ, KIT, ETH, MeteoSwiss  
!> Contact information: <https://icon-model.org>
!>
!> See AUTHORS.TXT for a list of authors  
!> See LICENSES/ for license information  
!> SPDX-License-Identifier: BSD-3-Clause
!> ------------------------------------------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### define assimilation config structure, read assimilation namelist and init configuration parameters
!>
MODULE mo_assimi_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, message_text, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_assimi_config

  ! ======================================================================================================= !
  !> ASSIMI_ configuration parameters
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_assimi_config
    LOGICAL                 :: init_running_means !< initialize npp buffer for NLCC process
    INTEGER                 :: ncanopy            !< number of canopy layers
    REAL(wp), ALLOCATABLE   :: canopy_bound_lai(:)
  CONTAINS
    PROCEDURE :: Init => Init_assimi_config
  END type t_assimi_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_config_class'

CONTAINS

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
    INTEGER                     :: jsbach_fixed_ncanopy       ! three canopy layer by default for MODEL_JSBACH
    INTEGER                     :: ncanopy_q_                 ! @TODO remove later, here and from namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    LOGICAL                     :: init_running_means
    INTEGER                     :: i
    INTEGER                     :: icanopy
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_canopy

    NAMELIST /jsb_assimi_nml/  &
      & active,                   &
      & ic_filename,              &
      & bc_filename,              &
      & ncanopy_q_,               &    ! @TODO remove later, here and from namelist
      & init_running_means

    INTEGER                     :: nml_handler, nml_unit, istat
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_assimi_config'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting assimi configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    jsbach_fixed_ncanopy          = 3            ! three canopy layer by default for MODEL_JSBACH
    ! Set namelist defaults
    active                        = .TRUE.
    bc_filename                   = 'bc_land_assimi.nc'
    ic_filename                   = 'ic_land_assimi.nc'
    init_running_means            = .FALSE.      ! .TRUE. .FALSE.
    ! ----------------------------------------------------------------------------------------------------- !
    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit    = position_nml('jsb_assimi_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_assimi_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Write namelist values into config
    config%active                       = active
    config%ic_filename                  = ic_filename
    config%bc_filename                  = bc_filename
    config%init_running_means           = init_running_means
    CALL close_nml(nml_handler)

    ! ncanopy is fixed to 3 layers
    ! for JSBACH ncanopy has so far not been a parameter in the assimi namelist but fixed to 3 in the code
    config%ncanopy = jsbach_fixed_ncanopy

    ! Create vertical axis for canopy layers
    vgrid_canopy  => new_vgrid('canopy_layer', ZAXIS_GENERIC, config%ncanopy, &
      & longname = 'Canopy layers from ASSIMI_ JSBACH', &
      & levels = (/ (REAL(i, kind=wp), i=1, config%ncanopy) /), units = '')
    CALL register_vgrid(vgrid_canopy)
    
    ! Prepare canopy boundaries of lai_cl
    ALLOCATE(config%canopy_bound_lai(0:config%ncanopy))
    DO icanopy = 0,config%ncanopy
      config%canopy_bound_lai(icanopy) = REAL(icanopy, wp) / REAL(config%ncanopy, wp)
      !write(*,*) "RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR: canopy_bound_lai(i)" , new%canopy_bound_lai(i)
    END DO

  END SUBROUTINE Init_assimi_config

#endif
END MODULE mo_assimi_config_class
