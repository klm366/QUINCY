!> Climate buffer process configuration
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>### Configuration of the climate buffer process
!>
!> The module defines the data structure [[t_climbuf_config]] containing the configuration of
!> process "CLIMBUF_", i.e. all climat buffer variables that are configurable via namelist.
!> In subroutine [[Init_climbuf_config]] we read the climate buffer namelist 'jsb_climbuf_nml'
!> and set the configuration data structure accordingly.
!>
MODULE mo_climbuf_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_climbuf_config

  !> The data structure contains all configurable parameters of process CLIMBUF. It extends
  !> [[mo_jsb_config_class:t_jsb_config]], a data structure containing basic parameters defined
  !> for all ICON-Land processes.
  !>
  !> _Note for the ford documentation_: The table lists the general parameters of
  !> [[mo_jsb_config_class:t_jsb_config]] as well as the parameters specific to the climate
  !> buffer actually defined in [[t_climbuf_config]].
  TYPE, EXTENDS(t_jsb_config) :: t_climbuf_config
    LOGICAL :: &
      & init_running_means                    !< Switch to initialize some running mean climate variables

   CONTAINS
     PROCEDURE :: Init => Init_climbuf_config
  END type t_climbuf_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> #### Setup the climate buffer configuration
  !
  !> In this subroutine we configure the climate buffer process based on the climate buffer namelist
  !> 'jsb_climbuf_nml'. Routine Init_climbuf_config is a type-bound procedure of [[t_climbuf_config]].
  !> It is called during initialization phase from [[mo_jsb_model_class:Configure_model_processes]]
  !> ('CALL this%processes(iproc)%p%Configure()').
  !>

  SUBROUTINE Init_climbuf_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_climbuf_config), INTENT(inout) :: config !<Configuration type of process (t_climbuf_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL ::                     &
      & active,                    & !< True: the CLIMBUF process is activated
      & lrestart_cont,             & !< Continue simulation although variables of the climbuf
                                     !< process are missing in restart file
      & init_running_means           !< Initialize running mean CLIMBUF variables

    CHARACTER(len=filename_max) :: &
      & ic_filename,               & !< Name of input file with initial conditions
      & bc_filename                  !< Name of input file with boundary conditions

    NAMELIST /jsb_climbuf_nml/     &
      & active,                    &
      & lrestart_cont,             &
      & ic_filename,               &
      & bc_filename,               &
      & init_running_means

    INTEGER ::                     &
      & nml_handler,               & !< Handle to deal with namelist file
      & nml_unit,                  & !< IO unit of the namelist file
      & istat                        !< Status flag

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_climbuf_config'
    ! -------------------------------------------------------------------------------------------------- !

    CALL message(TRIM(routine), 'Starting climbuf configuration')

    ! Set defaults
    active              = .FALSE.
    lrestart_cont       = .FALSE.      ! TRUE: Continue although climbuf variables are missing in restart file
    bc_filename         = 'bc_land_climbuf.nc'
    ic_filename         = 'ic_land_climbuf.nc'
    init_running_means  = .FALSE.

    ! Open the file containing the climate buffer namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    ! Find and read the climate buffer namelist
    nml_unit = position_nml('jsb_climbuf_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_climbuf_nml)

    CALL close_nml(nml_handler)

    !> After reading the namelist, all configurable climate buffer parameters - default parameters
    !> and values set in the namelist - are stored in the configuration data structure 'config' to
    !> make them easily accessible in all ICON-Land process modules ("dsl4jsb_Use_config(CLIMBUF_)").
    !>
    !> The specific choice of parameters is printed to the experiment's log file to document the
    !> setup.

    config%active              = active
    config%lrestart_cont       = lrestart_cont
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%init_running_means  = init_running_means

    IF (active .AND. init_running_means) CALL message(TRIM(routine), 'Initializing climate buffer')

  END SUBROUTINE Init_climbuf_config

#endif
END MODULE mo_climbuf_config_class
