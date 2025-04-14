!> Contains structures and methods for nix config
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

MODULE mo_nix_config_class

#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_nix_config

  ! ======================================================================================================= !
  !> NIX_ configuration parameters
  !>

  TYPE, EXTENDS(t_jsb_config) :: t_nix_config

!  ! Logicals
!  LOGICAL ::    &
!    active                    !< TRUE:  NIX is switched on 
                               !< FALSE: NIX is switched off     

  INTEGER ::    &
    nlev_sn                     !< Maximum number of snow layers



  CONTAINS

    PROCEDURE :: Init => Init_nix_config

  END type t_nix_config


  CHARACTER(len=*), PARAMETER :: modname = 'mo_nix_config_class'


CONTAINS

  ! ======================================================================================================= !
  !> read NIX_ namelist and init configuration parameters
  !>

  SUBROUTINE Init_nix_config(config)

    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_nix_config), INTENT(inout) :: config

    ! ----------------------------------------------------------------------------------------------------- !

    LOGICAL :: active

    INTEGER :: nlev_sn

    NAMELIST /jsb_nix_nml/ &
      & active,                &
      & nlev_sn                 


    INTEGER :: nml_handler, nml_unit, istat

    TYPE(t_jsb_vgrid), POINTER :: snow_vgrid

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_nix_config'

    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Starting NIX configuration')
    
      ! ----------------------------------------------------------------------------------------------------- !
      ! Set defaults
      active          = .FALSE.
      nlev_sn         = 10


      ! ----------------------------------------------------------------------------------------------------- !
      ! Read namelist
      nml_handler = open_nml(TRIM(config%namelist_filename))

      nml_unit = position_nml('jsb_nix_nml', nml_handler, STATUS=istat)
      IF (istat == POSITIONED) READ(nml_unit, jsb_nix_nml)

      CALL close_nml(nml_handler)
      
      ! ----------------------------------------------------------------------------------------------------- !
      ! Write namelist values into config FIXME: This would also be the place for call's to finish incase settings don't work
      !                                          compare sse_config
        config%active     = active
        config%nlev_sn    = nlev_sn 

      ! ----------------------------------------------------------------------------------------------------- !
      ! Add NIX specific vertical grid

      snow_vgrid  => new_vgrid('nix_vgrid', ZAXIS_GENERIC, nlev_sn, units='m')
      
      CALL register_vgrid(snow_vgrid)

       
  END SUBROUTINE Init_nix_config


#endif

END MODULE mo_nix_config_class

