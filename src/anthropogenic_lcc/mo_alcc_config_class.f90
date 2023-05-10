!> Summary: alcc config - defines vgrid required in alcc memory
!>
!>### Could contain namelist info (e.g. alcc scheme), but currently of no use beyond infrastructural requirements
!>
!>
!> Copyright (c) 2018 Rainer Schneck
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version   by Rainer Schneck (2018-11-28)
!> condensed       by Julia Nabel    (2021-06-23)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_alcc_config_class
#ifndef __NO_JSBACH__
  
  USE mo_exception,          ONLY: message, finish
  USE mo_io_units,           ONLY: filename_max
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_kind,               ONLY: wp
  USE mo_jsb_config_class,   ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_alcc_config

  TYPE, EXTENDS(t_jsb_config)      :: t_alcc_config
     INTEGER                       :: nr_of_pfts
     CHARACTER(len=filename_max-7) :: alcc_filename_prefix 
     CHARACTER(len=15)             :: scheme
     INTEGER                       :: interactive_istep
     INTEGER                       :: interactive_ipft_from
     INTEGER                       :: interactive_ipft_to
     REAL(wp)                      :: interactive_center(2)
     REAL(wp)                      :: interactive_radius
     REAL(wp)                      :: interactive_fract_max ! Maximum tile fraction on tile specified by interactive_ipft_from
     !< Prefix of alcc filename prefix with len filename_max-7 (to add yyyy.nc) 
   CONTAINS
     PROCEDURE :: Init => Init_alcc_config
  END type t_alcc_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize alcc process
  !
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize alcc process
  !!
  !! @param[inout]     config     Configuration type of process (t_alcc_config)
  !!  
  SUBROUTINE Init_alcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_alcc_config), INTENT(inout) :: config !<Configuration type of process (t_alcc_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                       :: active
    CHARACTER(len=filename_max)   :: ic_filename, bc_filename
    CHARACTER(len=filename_max-7) :: alcc_filename_prefix
    CHARACTER(len=15)             :: scheme
    TYPE(t_jsb_vgrid), POINTER    :: pft_layer
    INTEGER                       :: i, nr_of_pfts
    INTEGER                       :: interactive_istep
    INTEGER                       :: interactive_ipft_from
    INTEGER                       :: interactive_ipft_to
    REAL(wp)                      :: interactive_center(2)
    REAL(wp)                      :: interactive_radius
    REAL(wp)                      :: interactive_fract_max

    NAMELIST /jsb_alcc_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename,                 &
         alcc_filename_prefix,        &
         scheme,                      &
         nr_of_pfts,                  &
         interactive_istep,           &
         interactive_ipft_from,       &
         interactive_ipft_to,         &
         interactive_center,          &
         interactive_radius,          &
         interactive_fract_max

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_alcc_config'
    ! -------------------------------------------------------------------------------------------------- !

    CALL message(TRIM(routine), 'Starting alcc configuration')

    ! Set defaults
    active                  = .FALSE.
    bc_filename             = 'bc_land_alcc.nc'
    ic_filename             = 'ic_land_alcc.nc'
    alcc_filename_prefix    = 'bc_land_frac_11pfts_'
    nr_of_pfts              = 11
    scheme                  = 'maps'
    interactive_istep       = -1
    interactive_ipft_from   = 1                 ! Tropical broadleaf evergreen
    interactive_ipft_to     = 11                ! C3 crops
    interactive_center(1:2) = [-65._wp, -5._wp] ! Somewhere in Amazon
    interactive_radius      = 1.E6_wp           ! 1000 km
    interactive_fract_max   = 0.1_wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_alcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_alcc_nml)

    CALL close_nml(nml_handler)

    config%active                = active
    config%ic_filename           = ic_filename
    config%bc_filename           = bc_filename
    config%alcc_filename_prefix  = alcc_filename_prefix
    config%scheme                = TRIM(scheme)
    config%nr_of_pfts            = nr_of_pfts
    config%interactive_istep     = interactive_istep
    config%interactive_ipft_from = interactive_ipft_from
    config%interactive_ipft_to   = interactive_ipft_to
    config%interactive_center(:) = interactive_center(:)
    config%interactive_radius    = interactive_radius
    config%interactive_fract_max = interactive_fract_max

    SELECT CASE (TRIM(config%scheme))
    CASE ('maps')
      CALL message(TRIM(routine), 'Reading land use data from maps.')
    CASE ('interactive')
      CALL message(TRIM(routine), 'Interactively changing land use data.')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Unimplemented alcc scheme specified: '// TRIM(config%scheme))
    END SELECT

    ! Create vertical axis
    ! JN-TODO (RS...): currently I get nr_of_pfts from namelist -> better way? usecase?
    ! -> + this grid might be useful in other contexts, too?
    pft_layer  => new_vgrid('pfts', ZAXIS_GENERIC, nr_of_pfts, &
      & levels=(/ (REAL(i,kind=wp),i=1,nr_of_pfts) /), units='')
    CALL register_vgrid(pft_layer)

  END SUBROUTINE Init_alcc_config

#endif
END MODULE mo_alcc_config_class
