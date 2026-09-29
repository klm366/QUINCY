!>
!> Contains structures and methods for the topmodel for calculating innundation fractions
!>
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
!>
!> author
!>  Stiig Wilkenskjeld, MPI-M Hamburg
!>
!> Revision History
!> Adopted from JSBACH3 code by Thomas Kleinen, 2024-01-23
!>

MODULE mo_megan_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: r2s => real2string, i2s => int2string, l2s => logical2string
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_megan_config

  TYPE, EXTENDS(t_jsb_config) :: t_megan_config
    LOGICAL  :: dbg_process          ! Debug topmodel processes and tasks
    REAL(wp) :: scale_iso            ! Emission factor scaling, isoprene
    REAL(wp) :: scale_co             ! Emission factor scaling, carbon monoxide
    REAL(wp) :: scale_no             ! Emission factor scaling, nitric oxide
    REAL(wp) :: beta_co              ! Temperature adjustment factor, carbon monoxide
    REAL(wp) :: beta_no              ! Temperature adjustment factor, nitric oxide
    REAL(wp) :: relax_short          ! Relaxation time for short relaxation [s]
    REAL(wp) :: relax_long           ! Relaxation time for long relaxation [s]
  CONTAINS
    PROCEDURE :: Init => Init_megan_config
  END TYPE t_megan_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_megan_config_class'

CONTAINS

  SUBROUTINE Init_megan_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input

    CLASS(t_megan_config), INTENT(inout) :: config

    ! Namelist variables
    LOGICAL                     :: active
    LOGICAL                     :: dbg_process
    REAL(wp)                    :: scale_ISO
    REAL(wp)                    :: scale_CO
    REAL(wp)                    :: scale_NO
    REAL(wp)                    :: beta_CO
    REAL(wp)                    :: beta_NO
    REAL(wp)                    :: relax_short_s ! unit [s]
    REAL(wp)                    :: relax_long_s  ! unit [s]
    REAL(wp)                    :: relax_short   ! unit [d]
    REAL(wp)                    :: relax_long    ! unit [d]
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_megan_nml/                         &
      & active, dbg_process,                         &
      & scale_ISO, scale_CO, scale_NO,               &
      & beta_CO, beta_NO,                            &
      & relax_short_s, relax_long_s,                 &
      & relax_short, relax_long,                     &
      & ic_filename, bc_filename

    ! Locals
    INTEGER :: nml_handler, nml_unit, istat, i
    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_megan_config'

    CALL message(TRIM(routine), 'Starting megan module configuration')

    ! Set defaults
    active               = .FALSE.
    dbg_process          = .FALSE.
    ic_filename          = 'ic_land_soil.nc'
    bc_filename          = 'bc_land_soil.nc'
    !---emission factor scaling
    scale_ISO            = 0.589_wp             ! isoprene
    scale_CO             = 0.863_wp             ! carbon monoxide
    scale_NO             = 2.174_wp             ! nitric oxide
    !---temperature adjustment factor
    beta_CO              = 0.08_wp              ! carbon monoxide
    beta_NO              = 0.11_wp              ! nitric oxide
    ! Relaxation times (doubly defined)
    relax_short_s        = -1._wp               ! Impossible value to indicate unused
    relax_long_s         = -1._wp               ! Impossible value to indicate unused
    relax_short          =  1._wp               ! "Short" is by default 1 day
    relax_long           = 30._wp               ! "Long" is by default 30 days

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_megan_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_megan_nml)

    CALL close_nml(nml_handler)

    config%active       = active
    config%dbg_process  = dbg_process
    config%ic_filename  = ic_filename
    config%bc_filename  = bc_filename
    config%scale_ISO    = scale_ISO
    config%scale_CO     = scale_CO
    config%scale_NO     = scale_NO
    config%beta_CO      = beta_CO
    config%beta_NO      = beta_NO
    IF (relax_short_s > 0._wp) THEN
      config%relax_short = relax_short_s
    ELSE
      config%relax_short = relax_short * 86400._wp
    ENDIF
    IF (relax_long_s > 0._wp) THEN
      config%relax_long = relax_long_s
    ELSE
      config%relax_long = relax_long * 86400._wp
    ENDIF

    ! Tests and messsages
    IF (active) THEN
        CALL message(routine,'Megan module configured')
    ELSE
      CALL message(routine,'Megan module is not used')
    ENDIF

    !> -----------------------------------------------------------------------------------------------------------------------------
    IF (.NOT. active) RETURN

  END SUBROUTINE Init_megan_config

#endif
END MODULE mo_megan_config_class
