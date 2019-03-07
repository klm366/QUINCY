!>
!! @brief Contains structures and methods for forest-age config
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2017-11-20)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
MODULE mo_fage_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_util,              ONLY: int2string
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_fage_config

  TYPE, EXTENDS(t_jsb_config) :: t_fage_config
    INTEGER                   :: max_age                        ! Maximum age tracked  
    INTEGER                   :: ncohorts                       ! number of age cohorts
    CHARACTER(len=20)         :: age_cohort_scheme
    CHARACTER(len=10)         :: age_cohort_area_fracts
    CHARACTER(len=20)         :: disturb_distribution_scheme
  CONTAINS
    PROCEDURE             :: Init => Init_fage_config
  END type t_fage_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_config_class'

CONTAINS

  SUBROUTINE Init_fage_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC 

    CLASS(t_fage_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    INTEGER                     :: max_age

    NAMELIST /jsb_fage_nml/         &
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & max_age,                    &
      & ncohorts,                   &
      & age_cohort_scheme,          &
      & age_cohort_area_fracts,     &
      & disturb_distribution_scheme

    INTEGER :: nml_handler, nml_unit, istat
    INTEGER :: i, ncohorts
    CHARACTER(len=20)     :: age_cohort_scheme
    CHARACTER(len=10)     :: age_cohort_area_fracts
    CHARACTER(len=20)     :: disturb_distribution_scheme
    REAL(wp), ALLOCATABLE :: cohorts_grid_levels(:)
    REAL(wp), ALLOCATABLE :: fage_grid_levels(:)
    TYPE(t_jsb_vgrid), POINTER :: forest_age, age_cohorts_grid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_fage_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting fage configuration')

    ! Set defaults
    active            = .FALSE.
    bc_filename       = 'bc_land_fage.nc'  !TODO-JN do I need it at all?
    ic_filename       = 'ic_land_fage.nc'  !TODO-JN
    max_age           = 150
    ncohorts          = 16
    age_cohort_scheme = "static-equal" ! alternative: "static-shorterYoung"
    !TODO-JN: dynamic-biomass, dynamic-MaxLai, ...)
    age_cohort_area_fracts = "allFirst"     !TODO-JN: fromFile; allLast
    disturb_distribution_scheme = "equal"

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_fage_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_fage_nml)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%max_age             = max_age
    config%ncohorts            = ncohorts
    config%age_cohort_scheme   = TRIM(age_cohort_scheme)
    config%age_cohort_area_fracts = TRIM(age_cohort_area_fracts)
    config%disturb_distribution_scheme = TRIM(disturb_distribution_scheme)

    ALLOCATE(cohorts_grid_levels(ncohorts))
    ALLOCATE(fage_grid_levels(max_age))

    !Assertion: 
    IF ( max_age .LT. ncohorts) CALL finish(TRIM(routine), 'Assertion: number of age cohorts needs to be <= max age!')

!TODO-JN-COHORTS assert ncohorts < X 

    CALL close_nml(nml_handler)

    IF (.NOT. active) RETURN

    ! vertical grid with cohort dimension
    age_cohorts_grid  => new_vgrid('age_cohorts_grid', ZAXIS_GENERIC, ncohorts, &
      & levels  = (/ (REAL(i,kind=wp),i=1,ncohorts) /), &
      & lbounds = (/ (REAL(i,kind=wp),i=1,ncohorts) /), &
      & ubounds = (/ (REAL(i,kind=wp),i=1,ncohorts) /), &
      & units='')
    CALL register_vgrid(age_cohorts_grid)

    ! vertical grid with tracked age dimension 
    fage_grid_levels(1:max_age) = (/ (REAL(i,kind=wp),i=1,max_age) /)

    forest_age  => new_vgrid('forest_age', ZAXIS_GENERIC, max_age, &
      & levels  = fage_grid_levels (1:max_age) - 0.5_wp, &
      & lbounds = fage_grid_levels (1:max_age) - 1._wp,  &
      & ubounds = fage_grid_levels (1:max_age),          &
      & units='y')

!    forest_age  => new_vgrid('forest_age', ZAXIS_GENERIC, max_age, &
!      levels=(/ (REAL(i,kind=wp),i=1,max_age) /), units='')
    CALL register_vgrid(forest_age)

    CALL message(TRIM(routine), 'Tracked forest ages [y]: '//TRIM(int2string(max_age)))
    CALL message(TRIM(routine), 'Run with n age cohorts: '//TRIM(int2string(ncohorts)))


    DEALLOCATE(fage_grid_levels)

  END SUBROUTINE Init_fage_config

#endif
END MODULE mo_fage_config_class
