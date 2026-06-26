!> Contains structures and methods for alcc config
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
!>#### Configuration of the alcc process
!>
MODULE mo_alcc_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,          ONLY: message, finish
  USE mo_io_units,           ONLY: filename_max
  USE mo_util,               ONLY: tolower
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
  USE mo_kind,               ONLY: wp
  USE mo_jsb_config_class,   ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_alcc_config

  TYPE, EXTENDS(t_jsb_config)      :: t_alcc_config
     INTEGER                       :: nr_of_pfts
     CHARACTER(len=filename_max-7) :: alcc_filename_prefix, alcc_transitions_filename_prefix
       !< Prefix of alcc filename prefix with len filename_max-7 (to add yyyy.nc)
     CHARACTER(len=SHORT_NAME_LEN) :: scheme
     CHARACTER(len=15), ALLOCATABLE :: trans_names(:)
        !< specifies the transition types used in case of 'transition' scheme
     INTEGER, ALLOCATABLE          :: from_class(:), to_class(:), pft_class(:)
        !< from_class and to_class specify which class (natural, pasture, crop) loses/gains area through a transition type
        !< pft_class specifies which class each PFT belongs to
     LOGICAL                       :: l_daily_alcc
       !< Determines if land use change is to be applied on a daily (TRUE) or annual (FALSE) basis.
       ! Note: This also determines the year of the land use maps file to read in: in case of "annual" the file for
       ! the current year is read, while in case of "daily" the next year's file is read as target land cover map.
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
  SUBROUTINE Init_alcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_alcc_config), INTENT(inout) :: config !<Configuration type of process (t_alcc_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                       :: active, lrestart_cont
    CHARACTER(len=filename_max)   :: ic_filename, bc_filename
    CHARACTER(len=filename_max-7) :: alcc_filename_prefix, alcc_transitions_filename_prefix
    CHARACTER(len=SHORT_NAME_LEN) :: scheme
    LOGICAL                       :: l_daily_alcc
    TYPE(t_jsb_vgrid), POINTER    :: pft_layer
    TYPE(t_jsb_vgrid), POINTER    :: trans_layer
    INTEGER                       :: i, nr_of_pfts
    ! trans_names, from_class, to_class and pft_class are filled below depending on use case (jsbach_pfts, jsbach_12pfts, etc.)
    CHARACTER(len=15), ALLOCATABLE :: trans_names(:)  !< names of transition types
    INTEGER, ALLOCATABLE           :: from_class(:)   !< class losing area per transition (1: natural, 2: pasture, 3: c3crop/crop, 4: c4crop)
    INTEGER, ALLOCATABLE           :: to_class(:)     !< class gaining area per transition
    INTEGER, ALLOCATABLE           :: pft_class(:)    !< class index for each PFT

    NAMELIST /jsb_alcc_nml/      &
         active,                 &
         lrestart_cont,          &
         ic_filename,            &
         bc_filename,            &
         alcc_filename_prefix,   &
         alcc_transitions_filename_prefix, &
         scheme,                 &
         l_daily_alcc

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_alcc_config'
    ! -------------------------------------------------------------------------------------------------- !

    CALL message(TRIM(routine), 'Starting alcc configuration')

    ! Set defaults
    active               = .FALSE.
    lrestart_cont        = .FALSE.       ! TRUE: Continue although ALCC variables are missing in restart file
    bc_filename          = 'bc_land_alcc.nc'
    ic_filename          = 'ic_land_alcc.nc'
    alcc_filename_prefix = 'bc_land_frac_11pfts_'
    alcc_transitions_filename_prefix = 'LUH_transitions_'
    nr_of_pfts           = 11
    scheme               = 'maps'        ! options: 'maps', 'transition'
    l_daily_alcc         = .TRUE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_alcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_alcc_nml)

    CALL close_nml(nml_handler)

    config%active               = active
    config%lrestart_cont        = lrestart_cont
    config%ic_filename          = ic_filename
    config%bc_filename          = bc_filename
    config%alcc_filename_prefix = alcc_filename_prefix
    config%alcc_transitions_filename_prefix = alcc_transitions_filename_prefix
    config%scheme               = scheme
    config%l_daily_alcc         = l_daily_alcc
    config%nr_of_pfts           = nr_of_pfts

    SELECT CASE (TRIM(tolower(config%scheme)))
    CASE ('transition')

      ! Fill trans_names, from_class, to_class and pft_class based on the model use case
      SELECT CASE (TRIM(config%model_config%usecase))
      CASE ('jsbach_pfts')
        ALLOCATE(trans_names(6))
        trans_names = [character(len=15) :: &             !< 6 transition types for 'jsbach_pfts' (11 PFTs)
          'crop2pasture',   'crop2natural',   &
          'pasture2crop',   'pasture2natural', &
          'natural2crop',   'natural2pasture']
        ALLOCATE(from_class(6), to_class(6))              ! Indices of classes losing/gaining area through a transition type
        from_class = [3, 3, 2, 2, 1, 1]                   ! 1: natural, 2: pasture, 3: crop
        to_class   = [2, 1, 3, 1, 3, 2]
        ALLOCATE(pft_class(11))
        pft_class  = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3]    ! 8x natural, 2x pasture, 1x crop
      CASE ('jsbach_12pfts')
        ALLOCATE(trans_names(12))
        trans_names = [character(len=15) :: &             !< 12 transition types for 'jsbach_12pfts' (12 PFTs)
          'c3crop2c4crop',  'c3crop2pasture',  'c3crop2natural',  &
          'c4crop2c3crop',  'c4crop2pasture',  'c4crop2natural',  &
          'pasture2c3crop', 'pasture2c4crop',  'pasture2natural', &
          'natural2c3crop', 'natural2c4crop',  'natural2pasture']
        ALLOCATE(from_class(12), to_class(12))            ! Indices of classes losing/gaining area through a transition type
        from_class = [3, 3, 3, 4, 4, 4, 2, 2, 2, 1, 1, 1] ! 1: natural, 2: pasture, 3: c3crop, 4: c4crop
        to_class   = [4, 2, 1, 3, 2, 1, 3, 4, 1, 3, 4, 2]
        ALLOCATE(pft_class(12))
        pft_class  = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 4] ! 8x natural, 2x pasture, 1x c3crop, 1x c4crop
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Transition setup not defined for usecase: '// TRIM(config%model_config%usecase))
      END SELECT

      config%trans_names          = trans_names
      config%from_class           = from_class
      config%to_class             = to_class
      config%pft_class            = pft_class

    END SELECT

    IF (config%active) THEN
      SELECT CASE (TRIM(tolower(config%scheme)))
      CASE ('maps')
        CALL message(TRIM(routine), 'Reading land use data from maps.')
        config%scheme=tolower(config%scheme)

        IF (l_daily_alcc) THEN
          CALL message(TRIM(routine), '... land use change will be applied on a daily basis, reading next years land cover map' &
            &                       //' as target land cover to reach by the end of the year.')
        ELSE
          CALL message(TRIM(routine), '... land use change will be applied on an annual basis,' &
            &                       //' reading the land cover map of the current year to be applied in the first timestep.')
        END IF

      CASE ('transition')
        CALL message(TRIM(routine), 'Reading land use change from transitions data.')
        DO i = 1, SIZE(trans_names)
          CALL message(TRIM(routine), 'Transition type (' // TRIM(ADJUSTL(trans_names(i))) // ')')
        END DO
        config%scheme=tolower(config%scheme)

        IF(l_daily_alcc) THEN
          CALL message(TRIM(routine), '... land use change will be applied on a daily basis.')
        ELSE
          CALL message(TRIM(routine), '... land use change will be applied on an annual basis.')
        END IF

        trans_layer  => new_vgrid('transitions', ZAXIS_GENERIC, SIZE(trans_names), &
          & levels=(/ (REAL(i,kind=wp),i=1,SIZE(trans_names)) /), units='')
        CALL register_vgrid(trans_layer)

      CASE DEFAULT
        CALL finish(TRIM(routine), 'Not implemented alcc scheme specified: '// TRIM(config%scheme))
      END SELECT
    END IF

  END SUBROUTINE Init_alcc_config

#endif
END MODULE mo_alcc_config_class
