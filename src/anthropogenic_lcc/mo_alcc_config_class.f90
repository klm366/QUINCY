!> alcc config - defines vgrid required in alcc memory
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
!>#### Could contain namelist info (e.g. alcc scheme), but currently of no use beyond infrastructural requirements
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

  TYPE, EXTENDS(t_jsb_config)        :: t_alcc_config
     INTEGER                         :: nr_of_pfts
     CHARACTER(len=filename_max-7)   :: alcc_filename_prefix
     CHARACTER(len=filename_max-7)   :: alcc_transition_filename_prefix   !Mingyue
     CHARACTER(len=12), DIMENSION(9) :: transitions                       !Mingyue add transition types of the C,N,P
     INTEGER,           DIMENSION(9) :: LC_type2LC_class_mapping          !Mingyue add mapping table of 11 pfts to C, N, P , where 1->crop, 2-> pasture, 3->natual
       !< Prefix of alcc filename prefix with len filename_max-7 (to add yyyy.nc)
     CHARACTER(len=SHORT_NAME_LEN)   :: scheme
     LOGICAL                         :: l_daily_alcc
       !< determines if land use change is to be applied on a daily (TRUE) or annual (FALSE) basis
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
    CHARACTER(len=filename_max-7) :: alcc_filename_prefix
    CHARACTER(len=filename_max-7) :: alcc_transition_filename_prefix   !Mingyue
    CHARACTER(len=12), DIMENSION(9) :: transitions                     !Mingyue
    INTEGER,           DIMENSION(9) :: LC_type2LC_class_mapping        !Mingyue
    CHARACTER(len=SHORT_NAME_LEN) :: scheme
    LOGICAL                       :: l_daily_alcc
    TYPE(t_jsb_vgrid), POINTER    :: pft_layer
    TYPE(t_jsb_vgrid), POINTER    :: trans_layer                       !Mingyue
    TYPE(t_jsb_vgrid), POINTER    :: trans_pft_layer                   !Mingyue
    INTEGER                       :: i, nr_of_pfts

    NAMELIST /jsb_alcc_nml/      &
         active,                 &
         lrestart_cont,          &
         ic_filename,            &
         bc_filename,            &
         alcc_filename_prefix,   &
         alcc_transition_filename_prefix,    &       !Mingyue
         transitions,                        &       !Mingyue
         LC_type2LC_class_mapping,           &       !Mingyue
         scheme,                 &
         l_daily_alcc,           &
         nr_of_pfts

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
    alcc_transition_filename_prefix =   'LUH_scaled_transitions_T63_'                                        ! Mingyue Zhang: Set defaults for transition data nc file prefix 
    transitions                     =   (/ 'crop2crop', 'crop2pasture', 'crop2natural',          &           ! Mingyue Zhang: defeine transition types end   
                                           'pasture2crop', 'pasture2pasture', 'pasture2natural', &
                                           'natural2crop', 'natural2pasture', 'natural2natural'  &   /)      ! TODO: here the transitions type should be print out                                                                                                     
    LC_type2LC_class_mapping = (/ 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1 /)                                        ! pft 1-8 (Natural) ! pft 9 (Pasture), pft11 (Crops) ! Mingyue Zhang:define land_cover mapping table
    nr_of_pfts           = 11
    l_daily_alcc         = .TRUE.
    scheme               = 'maps'  ! switch to transitions when using transition scheme (Mingyue Zhang)

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_alcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_alcc_nml)

    CALL close_nml(nml_handler)

    config%active                                    = active
    config%lrestart_cont                             = lrestart_cont
    config%ic_filename                               = ic_filename
    config%bc_filename                               = bc_filename
    config%alcc_filename_prefix                      = alcc_filename_prefix          
    config%transitions                               = transitions                                   ! Mingyue Zhang: add config for transitions type 
    config%alcc_transition_filename_prefix           = alcc_transition_filename_prefix               ! Mingyue Zhang: add config for transition data nc file prefix 
    config%LC_type2LC_class_mapping                  = LC_type2LC_class_mapping                      ! Mingyue Zhang: add config for transition table
    config%scheme                                    = scheme
    config%l_daily_alcc                              = l_daily_alcc
    config%nr_of_pfts                                = nr_of_pfts

    SELECT CASE (TRIM(config%scheme))
    CASE ('maps')
      CALL message(TRIM(routine), 'Reading land use data from maps.')
    CASE ('transitions')
      CALL message(TRIM(routine), 'Reading land use change data from transition matrix.')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Not implemented alcc scheme specified: '// TRIM(config%scheme))
    END SELECT

    IF(l_daily_alcc) THEN
      CALL message(TRIM(routine), '... land use change will be applied on a daily basis.')
    ELSE
      CALL message(TRIM(routine), '... land use change will be applied on an annual basis.')
    END IF

    ! Create vertical axis
    ! JN-TODO (RS...): currently I get nr_of_pfts from namelist -> better way? usecase?
    ! -> + this grid might be useful in other contexts, too?
    pft_layer  => new_vgrid('pfts', ZAXIS_GENERIC, nr_of_pfts, &
      & levels=(/ (REAL(i,kind=wp),i=1,nr_of_pfts) /), units='')
    CALL register_vgrid(pft_layer)
    
    ! Mingyue Zhang define v_grid of (c,n,p) transition
    trans_layer  => new_vgrid('transition', ZAXIS_GENERIC, 9, &
      & levels=(/ (REAL(i,kind=wp),i=1,9) /), units='')
    CALL register_vgrid(trans_layer)

    ! Mingyue Zhang define v_grid of pfts transition
    trans_pft_layer  => new_vgrid('transition_pfts', ZAXIS_GENERIC, nr_of_pfts*nr_of_pfts, &
      & levels=(/ (REAL(i,kind=wp),i=1,nr_of_pfts*nr_of_pfts) /), units='')
    CALL register_vgrid(trans_pft_layer)

  END SUBROUTINE Init_alcc_config

#endif
END MODULE mo_alcc_config_class
