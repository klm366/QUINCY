!> QUINCY helper routines for model output as ascii text
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### various helper routines for the model output configuration
!>
MODULE mo_qs_output_txt_util
#ifndef __NO_QUINCY__

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message, message_text, finish

  IMPLICIT NONE
  PUBLIC        ! all is public

  ! ====================================================================================================== !
  !>
  !> write status and basename of txt output files
  !> status: 'is_active = TRUE' means that model output enabled for the particular file
  !>
  ! vegetation and soil pools
  CHARACTER(len=*), PARAMETER :: name_vegpool_c                         = 'vegpool_c'
  LOGICAL                     :: is_active_txt_vegpool_c                = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegpool_n                         = 'vegpool_n'
  LOGICAL                     :: is_active_txt_vegpool_n                = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegpool_p                         = 'vegpool_p'
  LOGICAL                     :: is_active_txt_vegpool_p                = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegpool_dc13dc14                  = 'vegpool_dc13dc14'
  LOGICAL                     :: is_active_txt_vegpool_dc13dc14         = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegpool_dn15                      = 'vegpool_dn15'
  LOGICAL                     :: is_active_txt_vegpool_dn15             = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_org                      = 'soilpool_org'
  LOGICAL                     :: is_active_txt_soilpool_org             = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_inorg                    = 'soilpool_inorg'
  LOGICAL                     :: is_active_txt_soilpool_inorg           = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_dc13dc14dn15             = 'soilpool_dc13dc14dn15'
  LOGICAL                     :: is_active_txt_soilpool_dc13dc14dn15    = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_c_org_sl                 = 'soilpool_c_org_sl'
  LOGICAL                     :: is_active_txt_soilpool_c_org_sl        = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_n_org_sl                 = 'soilpool_n_org_sl'
  LOGICAL                     :: is_active_txt_soilpool_n_org_sl        = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_p_org_sl                 = 'soilpool_p_org_sl'
  LOGICAL                     :: is_active_txt_soilpool_p_org_sl        = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_ndn15_inorg_sl           = 'soilpool_ndn15_inorg_sl'
  LOGICAL                     :: is_active_txt_soilpool_ndn15_inorg_sl  = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_p_inorg_sl               = 'soilpool_p_inorg_sl'
  LOGICAL                     :: is_active_txt_soilpool_p_inorg_sl      = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_dc13_org_sl              = 'soilpool_dc13_org_sl'
  LOGICAL                     :: is_active_txt_soilpool_dc13_org_sl     = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_dc14_org_sl              = 'soilpool_dc14_org_sl'
  LOGICAL                     :: is_active_txt_soilpool_dc14_org_sl     = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilpool_dn15_org_sl              = 'soilpool_dn15_org_sl'
  LOGICAL                     :: is_active_txt_soilpool_dn15_org_sl     = .FALSE.
  ! ecosystem pool
  CHARACTER(len=*), PARAMETER :: name_ecosyspool_cnp          = 'ecosyspool_cnp'
  LOGICAL                     :: is_active_txt_ecosyspool_cnp = .FALSE.
  ! vegetation and soil fluxes
  CHARACTER(len=*), PARAMETER :: name_vegflux_c                   = 'vegflux_c'
  LOGICAL                     :: is_active_txt_vegflux_c          = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegflux_n                   = 'vegflux_n'
  LOGICAL                     :: is_active_txt_vegflux_n          = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegflux_p                   = 'vegflux_p'
  LOGICAL                     :: is_active_txt_vegflux_p          = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegflux_dc13dc14            = 'vegflux_dc13dc14'
  LOGICAL                     :: is_active_txt_vegflux_dc13dc14   = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegflux_h2o                 = 'vegflux_h2o'
  LOGICAL                     :: is_active_txt_vegflux_h2o        = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_vegflux_energy              = 'vegflux_energy'
  LOGICAL                     :: is_active_txt_vegflux_energy     = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilflux_cdc13dc14          = 'soilflux_cdc13dc14'
  LOGICAL                     :: is_active_txt_soilflux_cdc13dc14 = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilflux_n                  = 'soilflux_n'
  LOGICAL                     :: is_active_txt_soilflux_n         = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilflux_p                  = 'soilflux_p'
  LOGICAL                     :: is_active_txt_soilflux_p         = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilflux_dn15               = 'soilflux_dn15'
  LOGICAL                     :: is_active_txt_soilflux_dn15      = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilflux_c_sl               = 'soilflux_c_sl'
  LOGICAL                     :: is_active_txt_soilflux_c_sl      = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilflux_n_sl               = 'soilflux_n_sl'
  LOGICAL                     :: is_active_txt_soilflux_n_sl      = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soilflux_p_sl               = 'soilflux_p_sl'
  LOGICAL                     :: is_active_txt_soilflux_p_sl      = .FALSE.
  ! mycorrhiza fluxes
  CHARACTER(len=*), PARAMETER :: name_mycflux_cnp           = 'mycflux_cnp'
  LOGICAL                     :: is_active_txt_mycflux_cnp  = .FALSE.
  ! soil physics
  CHARACTER(len=*), PARAMETER :: name_soil_physics              = 'soil_physics'
  LOGICAL                     :: is_active_txt_soil_physics     = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soil_physics_sl           = 'soil_physics_sl'
  LOGICAL                     :: is_active_txt_soil_physics_sl  = .FALSE.
  ! diagnostics
  CHARACTER(len=*), PARAMETER :: name_veg_diagnostics                         = 'veg_diagnostics'
  LOGICAL                     :: is_active_txt_veg_diagnostics                = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_veg_canopy_diagnostics_cl               = 'veg_canopy_diagnostics_cl'
  LOGICAL                     :: is_active_txt_veg_canopy_diagnostics_cl      = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_myc_diagnostics                         = 'myc_diagnostics'
  LOGICAL                     :: is_active_txt_myc_diagnostics                = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soil_diagnostics_microbial_sl           = 'soil_diagnostics_microbial_sl'
  LOGICAL                     :: is_active_txt_soil_diagnostics_microbial_sl  = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_soil_diagnostics_sorption_sl            = 'soil_diagnostics_sorption_sl'
  LOGICAL                     :: is_active_txt_soil_diagnostics_sorption_sl   = .FALSE.
  ! matter conservation / mass balance
  CHARACTER(len=*), PARAMETER :: name_matter_conservation_test          = 'matter_conservation_test'
  LOGICAL                     :: is_active_txt_matter_conservation_test = .FALSE.
  CHARACTER(len=*), PARAMETER :: name_sb_mass_balance                   = 'sb_mass_balance'
  LOGICAL                     :: is_active_txt_sb_mass_balance          = .FALSE.
  ! atmospheric forcing
  CHARACTER(len=*), PARAMETER :: name_atm_forcing           = 'atm_forcing'
  LOGICAL                     :: is_active_txt_atm_forcing  = .FALSE.

  ! ====================================================================================================== !
  !>
  !> LUN of txt output files
  !>
  INTEGER, PARAMETER :: first_lun_ascii_txt_output_files = 100
  ENUM, BIND(C)
  ENUMERATOR :: &
    & lun_vegpool_c = first_lun_ascii_txt_output_files, &
    & lun_vegpool_n, &
    & lun_vegpool_p, &
    & lun_vegpool_dc13dc14, &
    & lun_vegpool_dn15, &
    & lun_soilpool_org, &
    & lun_soilpool_inorg, &
    & lun_soilpool_dc13dc14dn15, &
    & lun_soilpool_c_org_sl, &
    & lun_soilpool_n_org_sl, &
    & lun_soilpool_p_org_sl, &
    & lun_soilpool_ndn15_inorg_sl, &
    & lun_soilpool_p_inorg_sl, &
    & lun_soilpool_dc13_org_sl, &
    & lun_soilpool_dc14_org_sl, &
    & lun_soilpool_dn15_org_sl, &
    & lun_ecosyspool_cnp, &
    & lun_vegflux_c, &
    & lun_vegflux_n, &
    & lun_vegflux_p, &
    & lun_vegflux_dc13dc14, &
    & lun_vegflux_h2o, &
    & lun_vegflux_energy, &
    & lun_soilflux_cdc13dc14, &
    & lun_soilflux_n, &
    & lun_soilflux_p, &
    & lun_soilflux_dn15, &
    & lun_soilflux_c_sl, &
    & lun_soilflux_n_sl, &
    & lun_soilflux_p_sl, &
    & lun_mycflux_cnp, &
    & lun_soil_physics, &
    & lun_soil_physics_sl, &
    & lun_veg_diagnostics, &
    & lun_veg_canopy_diagnostics_cl, &
    & lun_myc_diagnostics, &
    & lun_soil_diagnostics_microbial_sl, &
    & lun_soil_diagnostics_sorption_sl, &
    & lun_matter_conservation_test, &
    & lun_sb_mass_balance, &
    & lun_atm_forcing, &
    & last_lun_ascii_txt_output_files_plus_one
  END ENUM

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_qs_output_txt_util'

CONTAINS

  ! ====================================================================================================== !
  !>
  !> print to logfile the configuration of ascii txt output files
  !>
  SUBROUTINE print_output_txt_config()

    CHARACTER(len=17),  PARAMETER :: format_message = '(A, I4, L2, A, A)'
    CHARACTER(len=*),   PARAMETER :: file_details   = '  LUN, active, filename: '
    CHARACTER(len=*),   PARAMETER :: routine = TRIM(modname)//':print_output_txt_config'

    WRITE(message_text, *) 'Show configuration of ascii txt output files:', NEW_LINE('a')
    CALL message('', message_text, all_print = .TRUE.)

    ! vegetation and soil pools
    WRITE(message_text, format_message) file_details, &
      &   lun_vegpool_c, is_active_txt_vegpool_c, ' ', name_vegpool_c
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text, format_message) file_details, &
      &   lun_vegpool_n, is_active_txt_vegpool_n, ' ', name_vegpool_n
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text, format_message) file_details, &
      &   lun_vegpool_p, is_active_txt_vegpool_p, ' ', name_vegpool_p
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_vegpool_dc13dc14, is_active_txt_vegpool_dc13dc14, ' ', name_vegpool_dc13dc14
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_vegpool_dn15, is_active_txt_vegpool_dn15, ' ', name_vegpool_dn15
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_org, is_active_txt_soilpool_org, ' ', name_soilpool_org
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_inorg, is_active_txt_soilpool_inorg, ' ', name_soilpool_inorg
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_dc13dc14dn15, is_active_txt_soilpool_dc13dc14dn15, ' ', name_soilpool_dc13dc14dn15
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_c_org_sl, is_active_txt_soilpool_c_org_sl, ' ', name_soilpool_c_org_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_n_org_sl, is_active_txt_soilpool_n_org_sl, ' ', name_soilpool_n_org_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_p_org_sl, is_active_txt_soilpool_p_org_sl, ' ', name_soilpool_p_org_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_ndn15_inorg_sl, is_active_txt_soilpool_ndn15_inorg_sl, ' ', name_soilpool_ndn15_inorg_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_p_inorg_sl, is_active_txt_soilpool_p_inorg_sl, ' ', name_soilpool_p_inorg_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_dc13_org_sl, is_active_txt_soilpool_dc13_org_sl, ' ', name_soilpool_dc13_org_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_dc14_org_sl, is_active_txt_soilpool_dc14_org_sl, ' ', name_soilpool_dc14_org_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilpool_dn15_org_sl, is_active_txt_soilpool_dn15_org_sl, ' ', name_soilpool_dn15_org_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    ! ecosystem pool
    WRITE(message_text,  format_message) file_details, &
      &   lun_ecosyspool_cnp, is_active_txt_ecosyspool_cnp, ' ', name_ecosyspool_cnp
    CALL message(' ', message_text, all_print = .TRUE.)
    ! vegetation and soil fluxes
    WRITE(message_text,  format_message) file_details, &
      &   lun_vegflux_c, is_active_txt_vegflux_c, ' ', name_vegflux_c
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_vegflux_n, is_active_txt_vegflux_n, ' ', name_vegflux_n
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_vegflux_p, is_active_txt_vegflux_p, ' ', name_vegflux_p
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_vegflux_dc13dc14, is_active_txt_vegflux_dc13dc14, ' ', name_vegflux_dc13dc14
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_vegflux_h2o, is_active_txt_vegflux_h2o, ' ', name_vegflux_h2o
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_vegflux_energy, is_active_txt_vegflux_energy, ' ', name_vegflux_energy
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilflux_cdc13dc14, is_active_txt_soilflux_cdc13dc14, ' ', name_soilflux_cdc13dc14
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilflux_n, is_active_txt_soilflux_n, ' ', name_soilflux_n
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilflux_p, is_active_txt_soilflux_p, ' ', name_soilflux_p
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilflux_dn15, is_active_txt_soilflux_dn15, ' ', name_soilflux_dn15
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilflux_c_sl, is_active_txt_soilflux_c_sl, ' ', name_soilflux_c_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilflux_n_sl, is_active_txt_soilflux_n_sl, ' ', name_soilflux_n_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soilflux_p_sl, is_active_txt_soilflux_p_sl, ' ', name_soilflux_p_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    ! mycorrhiza fluxes
    WRITE(message_text,  format_message) file_details, &
      &   lun_mycflux_cnp, is_active_txt_mycflux_cnp, ' ', name_mycflux_cnp
    CALL message(' ', message_text, all_print = .TRUE.)
    ! soil physics
    WRITE(message_text,  format_message) file_details, &
      &   lun_soil_physics, is_active_txt_soil_physics, ' ', name_soil_physics
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soil_physics_sl, is_active_txt_soil_physics_sl, ' ', name_soil_physics_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    ! diagnostics
    WRITE(message_text,  format_message) file_details, &
      &   lun_veg_diagnostics, is_active_txt_veg_diagnostics, ' ', name_veg_diagnostics
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_veg_canopy_diagnostics_cl, is_active_txt_veg_canopy_diagnostics_cl, ' ', name_veg_canopy_diagnostics_cl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_myc_diagnostics, is_active_txt_myc_diagnostics, ' ', name_myc_diagnostics
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soil_diagnostics_microbial_sl, is_active_txt_soil_diagnostics_microbial_sl, ' ', name_soil_diagnostics_microbial_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_soil_diagnostics_sorption_sl, is_active_txt_soil_diagnostics_sorption_sl, ' ', name_soil_diagnostics_sorption_sl
    CALL message(' ', message_text, all_print = .TRUE.)
    ! matter conservation / mass balance
    WRITE(message_text,  format_message) file_details, &
      &   lun_matter_conservation_test, is_active_txt_matter_conservation_test, ' ', name_matter_conservation_test
    CALL message(' ', message_text, all_print = .TRUE.)
    WRITE(message_text,  format_message) file_details, &
      &   lun_sb_mass_balance, is_active_txt_sb_mass_balance, ' ', name_sb_mass_balance
    CALL message(' ', message_text, all_print = .TRUE.)
    ! atmospheric forcing
    WRITE(message_text,  format_message) file_details, &
      &   lun_atm_forcing, is_active_txt_atm_forcing, ' ', name_atm_forcing
    CALL message(' ', message_text, all_print = .TRUE.)

    WRITE(message_text, *) TRIM(routine), ' Finished', NEW_LINE('a')
    CALL message('', message_text, all_print = .TRUE.)
  END SUBROUTINE print_output_txt_config

  ! ====================================================================================================== !
  !>
  !>configure the write-status of ascii txt output files
  !>
  !>  lun of files that are an argument of init_output_txt_file_set() are activated for output
  !>
  SUBROUTINE init_output_txt_config(output_txt_file_set)

    CHARACTER(*), INTENT(IN) :: output_txt_file_set

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_output_txt_config'

    ! write to logfile
    CALL message('', NEW_LINE('a'), all_print = .TRUE.)
    WRITE(message_text, *) TRIM(routine), ' Applying ascii txt output fileset: ', output_txt_file_set
    CALL message('', message_text, all_print = .TRUE.)
    ! write to screen/terminal
    PRINT *, NEW_LINE('a'), TRIM(message_text)

    ! init output fileset based on namelist option
    SELECT CASE (TRIM(output_txt_file_set))
    ! all: land refruns
    CASE ('all')
      CALL init_output_txt_file_set_all()

    ! basic: land
    CASE ('basic')
      CALL init_output_txt_file_set( &
        &   lun_vegpool_c, &
        &   lun_vegpool_n, &
        &   lun_vegpool_p, &
        &   lun_soilpool_org, &
        &   lun_soilpool_inorg, &
        &   lun_ecosyspool_cnp, &
        &   lun_vegflux_c, &
        &   lun_vegflux_n, &
        &   lun_vegflux_p, &
        &   lun_vegflux_h2o, &
        &   lun_vegflux_energy, &
        &   lun_soilflux_cdc13dc14, &
        &   lun_soilflux_n, &
        &   lun_soilflux_p, &
        &   lun_mycflux_cnp, &
        &   lun_soil_physics, &
        &   lun_soil_physics_sl, &
        &   lun_veg_diagnostics, &
        &   lun_veg_canopy_diagnostics_cl, &
        &   lun_matter_conservation_test, &
        &   lun_sb_mass_balance &
        &   )

    ! minimal: land, plant
    CASE ('minimal')
      CALL init_output_txt_file_set( &
        &   lun_vegpool_c, &
        &   lun_vegpool_n, &
        &   lun_vegpool_p, &
        &   lun_soilpool_org, &
        &   lun_soilpool_inorg, &
        &   lun_ecosyspool_cnp, &
        &   lun_vegflux_c, &
        &   lun_vegflux_h2o, &
        &   lun_vegflux_energy, &
        &   lun_veg_diagnostics, &
        &   lun_matter_conservation_test, &
        &   lun_sb_mass_balance &
        &   )

    ! plant
    CASE ('plant')
      CALL init_output_txt_file_set( &
        &   lun_vegpool_c, &
        &   lun_vegpool_n, &
        &   lun_vegpool_p, &
        &   lun_vegpool_dc13dc14, &
        &   lun_vegpool_dn15, &
        &   lun_soilpool_org, &
        &   lun_soilpool_inorg, &
        &   lun_soilpool_dc13dc14dn15, &
        &   lun_ecosyspool_cnp, &
        &   lun_vegflux_c, &
        &   lun_vegflux_n, &
        &   lun_vegflux_p, &
        &   lun_vegflux_dc13dc14, &
        &   lun_vegflux_h2o, &
        &   lun_vegflux_energy, &
        &   lun_soil_physics, &
        &   lun_soil_physics_sl, &
        &   lun_veg_diagnostics, &
        &   lun_veg_canopy_diagnostics_cl, &
        &   lun_atm_forcing &
        &   )

    ! canopy
    CASE ('canopy')
      CALL init_output_txt_file_set( &
        &   lun_vegflux_c, &
        &   lun_vegflux_dc13dc14, &
        &   lun_vegflux_h2o, &
        &   lun_vegflux_energy, &
        &   lun_veg_canopy_diagnostics_cl, &
        &   lun_atm_forcing &
        &   )

    ! soil
    CASE ('soil')
      CALL init_output_txt_file_set( &
        &   lun_soilpool_org, &
        &   lun_soilpool_c_org_sl, &
        &   lun_soilpool_n_org_sl, &
        &   lun_soilpool_p_org_sl, &
        &   lun_soilpool_ndn15_inorg_sl, &
        &   lun_soilpool_p_inorg_sl, &
        &   lun_soilpool_dc13_org_sl, &
        &   lun_soilpool_dc14_org_sl, &
        &   lun_soilpool_dn15_org_sl, &
        &   lun_soilflux_cdc13dc14, &
        &   lun_soilflux_n, &
        &   lun_soilflux_p, &
        &   lun_soilflux_dn15, &
        &   lun_soilflux_c_sl, &
        &   lun_soilflux_n_sl, &
        &   lun_soilflux_p_sl, &
        &   lun_soil_diagnostics_microbial_sl, &
        &   lun_soil_diagnostics_sorption_sl, &
        &   lun_sb_mass_balance &
        &   )

    ! catch wrong setting for output fileset
    CASE DEFAULT
      CALL finish(TRIM(routine), 'No such ascii txt output fileset defined.')
    END SELECT

    WRITE(message_text, *) TRIM(routine), ' Finished', NEW_LINE('a')
    CALL message('', message_text, all_print = .TRUE.)
  END SUBROUTINE init_output_txt_config

  ! ====================================================================================================== !
  !>
  !> init all ascii txt output files
  !>
  SUBROUTINE init_output_txt_file_set_all()

    INTEGER                     :: ilun
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_output_txt_file_set_all'

    DO ilun = first_lun_ascii_txt_output_files, (last_lun_ascii_txt_output_files_plus_one - 1)
      CALL set_active_output_txt_file(ilun)
    END DO
  END SUBROUTINE init_output_txt_file_set_all

  ! ====================================================================================================== !
  !>
  !> init a specific set of ascii txt output files
  !>
  SUBROUTINE init_output_txt_file_set( &
    &   lun1,  lun2,  lun3,  lun4,  lun5,  lun6,  lun7,  lun8,  lun9,  lun10, &
    &   lun11, lun12, lun13, lun14, lun15, lun16, lun17, lun18, lun19, lun20, &
    &   lun21, lun22, lun23, lun24, lun25, lun26, lun27, lun28, lun29, lun30, &
    &   lun31, lun32, lun33, lun34, lun35, lun36, lun37, lun38, lun39, lun40, &
    &   lun41, lun42, lun43, lun44, lun45, lun46, lun47, lun48, lun49, lun50 &
    &   )

    INTEGER, INTENT(IN), OPTIONAL :: &
      &   lun1,  lun2,  lun3,  lun4,  lun5,  lun6,  lun7,  lun8,  lun9,  lun10, &
      &   lun11, lun12, lun13, lun14, lun15, lun16, lun17, lun18, lun19, lun20, &
      &   lun21, lun22, lun23, lun24, lun25, lun26, lun27, lun28, lun29, lun30, &
      &   lun31, lun32, lun33, lun34, lun35, lun36, lun37, lun38, lun39, lun40, &
      &   lun41, lun42, lun43, lun44, lun45, lun46, lun47, lun48, lun49, lun50

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_output_txt_file_set'

    IF (PRESENT (lun1))  CALL set_active_output_txt_file(lun1)
    IF (PRESENT (lun2))  CALL set_active_output_txt_file(lun2)
    IF (PRESENT (lun3))  CALL set_active_output_txt_file(lun3)
    IF (PRESENT (lun4))  CALL set_active_output_txt_file(lun4)
    IF (PRESENT (lun5))  CALL set_active_output_txt_file(lun5)
    IF (PRESENT (lun6))  CALL set_active_output_txt_file(lun6)
    IF (PRESENT (lun7))  CALL set_active_output_txt_file(lun7)
    IF (PRESENT (lun8))  CALL set_active_output_txt_file(lun8)
    IF (PRESENT (lun9))  CALL set_active_output_txt_file(lun9)
    IF (PRESENT (lun10)) CALL set_active_output_txt_file(lun10)
    IF (PRESENT (lun11)) CALL set_active_output_txt_file(lun11)
    IF (PRESENT (lun12)) CALL set_active_output_txt_file(lun12)
    IF (PRESENT (lun13)) CALL set_active_output_txt_file(lun13)
    IF (PRESENT (lun14)) CALL set_active_output_txt_file(lun14)
    IF (PRESENT (lun15)) CALL set_active_output_txt_file(lun15)
    IF (PRESENT (lun16)) CALL set_active_output_txt_file(lun16)
    IF (PRESENT (lun17)) CALL set_active_output_txt_file(lun17)
    IF (PRESENT (lun18)) CALL set_active_output_txt_file(lun18)
    IF (PRESENT (lun19)) CALL set_active_output_txt_file(lun19)
    IF (PRESENT (lun20)) CALL set_active_output_txt_file(lun20)
    IF (PRESENT (lun21)) CALL set_active_output_txt_file(lun21)
    IF (PRESENT (lun22)) CALL set_active_output_txt_file(lun22)
    IF (PRESENT (lun23)) CALL set_active_output_txt_file(lun23)
    IF (PRESENT (lun24)) CALL set_active_output_txt_file(lun24)
    IF (PRESENT (lun25)) CALL set_active_output_txt_file(lun25)
    IF (PRESENT (lun26)) CALL set_active_output_txt_file(lun26)
    IF (PRESENT (lun27)) CALL set_active_output_txt_file(lun27)
    IF (PRESENT (lun28)) CALL set_active_output_txt_file(lun28)
    IF (PRESENT (lun29)) CALL set_active_output_txt_file(lun29)
    IF (PRESENT (lun30)) CALL set_active_output_txt_file(lun30)
    IF (PRESENT (lun31)) CALL set_active_output_txt_file(lun31)
    IF (PRESENT (lun32)) CALL set_active_output_txt_file(lun32)
    IF (PRESENT (lun33)) CALL set_active_output_txt_file(lun33)
    IF (PRESENT (lun34)) CALL set_active_output_txt_file(lun34)
    IF (PRESENT (lun35)) CALL set_active_output_txt_file(lun35)
    IF (PRESENT (lun36)) CALL set_active_output_txt_file(lun36)
    IF (PRESENT (lun37)) CALL set_active_output_txt_file(lun37)
    IF (PRESENT (lun38)) CALL set_active_output_txt_file(lun38)
    IF (PRESENT (lun39)) CALL set_active_output_txt_file(lun39)
    IF (PRESENT (lun40)) CALL set_active_output_txt_file(lun40)
    IF (PRESENT (lun41)) CALL set_active_output_txt_file(lun41)
    IF (PRESENT (lun42)) CALL set_active_output_txt_file(lun42)
    IF (PRESENT (lun43)) CALL set_active_output_txt_file(lun43)
    IF (PRESENT (lun44)) CALL set_active_output_txt_file(lun44)
    IF (PRESENT (lun45)) CALL set_active_output_txt_file(lun45)
    IF (PRESENT (lun46)) CALL set_active_output_txt_file(lun46)
    IF (PRESENT (lun47)) CALL set_active_output_txt_file(lun47)
    IF (PRESENT (lun48)) CALL set_active_output_txt_file(lun48)
    IF (PRESENT (lun49)) CALL set_active_output_txt_file(lun49)
    IF (PRESENT (lun50)) CALL set_active_output_txt_file(lun50)
  END SUBROUTINE init_output_txt_file_set

  ! ====================================================================================================== !
  !>
  !> activate output for a specific ascii txt output file
  !>
  !> test if an output file is matching the LUN
  !>
  SUBROUTINE set_active_output_txt_file(lun_test)

    INTEGER, INTENT(IN) :: lun_test     !< LUN of the output file to activate output writing

    LOGICAL                     :: is_valid_lun
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_active_output_txt_file'

    ! catch LUN outside of defined values
    IF (lun_test < first_lun_ascii_txt_output_files .OR. lun_test > (last_lun_ascii_txt_output_files_plus_one - 1)) THEN
      WRITE(message_text, '(A,I4,A,I4,A,I5)') ' The LUN you defined for the output txt file is outside of defined values (', &
        & first_lun_ascii_txt_output_files, ',', (last_lun_ascii_txt_output_files_plus_one - 1), '):', lun_test
      CALL finish(TRIM(routine), message_text)
    END IF

    ! set test variable to FALSE
    is_valid_lun = .FALSE.

    ! activate the specific ascii txt output file that has the LUN 'lun_test'
    IF (lun_test == lun_vegpool_c) THEN
      is_active_txt_vegpool_c = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegpool_n) THEN
      is_active_txt_vegpool_n = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegpool_p) THEN
      is_active_txt_vegpool_p = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegpool_dc13dc14) THEN
      is_active_txt_vegpool_dc13dc14 = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegpool_dn15) THEN
      is_active_txt_vegpool_dn15 = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_org) THEN
      is_active_txt_soilpool_org = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_inorg) THEN
      is_active_txt_soilpool_inorg = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_dc13dc14dn15) THEN
      is_active_txt_soilpool_dc13dc14dn15 = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_c_org_sl) THEN
      is_active_txt_soilpool_c_org_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_n_org_sl) THEN
      is_active_txt_soilpool_n_org_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_p_org_sl) THEN
      is_active_txt_soilpool_p_org_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_ndn15_inorg_sl) THEN
      is_active_txt_soilpool_ndn15_inorg_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_p_inorg_sl) THEN
      is_active_txt_soilpool_p_inorg_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_dc13_org_sl) THEN
      is_active_txt_soilpool_dc13_org_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_dc14_org_sl) THEN
      is_active_txt_soilpool_dc14_org_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilpool_dn15_org_sl) THEN
      is_active_txt_soilpool_dn15_org_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_ecosyspool_cnp) THEN
      is_active_txt_ecosyspool_cnp = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegflux_c) THEN
      is_active_txt_vegflux_c = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegflux_n) THEN
      is_active_txt_vegflux_n = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegflux_p) THEN
      is_active_txt_vegflux_p = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegflux_dc13dc14) THEN
      is_active_txt_vegflux_dc13dc14 = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegflux_h2o) THEN
      is_active_txt_vegflux_h2o = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_vegflux_energy) THEN
      is_active_txt_vegflux_energy = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilflux_cdc13dc14) THEN
      is_active_txt_soilflux_cdc13dc14 = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilflux_n) THEN
      is_active_txt_soilflux_n = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilflux_p) THEN
      is_active_txt_soilflux_p = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilflux_dn15) THEN
      is_active_txt_soilflux_dn15 = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilflux_c_sl) THEN
      is_active_txt_soilflux_c_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilflux_n_sl) THEN
      is_active_txt_soilflux_n_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soilflux_p_sl) THEN
      is_active_txt_soilflux_p_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_mycflux_cnp) THEN
      is_active_txt_mycflux_cnp = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soil_physics) THEN
      is_active_txt_soil_physics = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soil_physics_sl) THEN
      is_active_txt_soil_physics_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_veg_diagnostics) THEN
      is_active_txt_veg_diagnostics = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_veg_canopy_diagnostics_cl) THEN
      is_active_txt_veg_canopy_diagnostics_cl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_myc_diagnostics) THEN
      is_active_txt_myc_diagnostics = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soil_diagnostics_microbial_sl) THEN
      is_active_txt_soil_diagnostics_microbial_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_soil_diagnostics_sorption_sl) THEN
      is_active_txt_soil_diagnostics_sorption_sl = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_matter_conservation_test) THEN
      is_active_txt_matter_conservation_test = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_sb_mass_balance) THEN
      is_active_txt_sb_mass_balance = .TRUE.
      is_valid_lun = .TRUE.
    END IF
    IF (lun_test == lun_atm_forcing) THEN
      is_active_txt_atm_forcing = .TRUE.
      is_valid_lun = .TRUE.
    END IF

    ! finish if no output file is matching the LUN
    IF (.NOT. is_valid_lun) THEN
      WRITE(message_text, *) ' No output file is matching the LUN: ', lun_test
      CALL finish(TRIM(routine), message_text)
    END IF
  END SUBROUTINE set_active_output_txt_file

#endif
END MODULE mo_qs_output_txt_util
