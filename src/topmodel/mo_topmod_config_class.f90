!>
!! @brief Contains structures and methods for the topmodel for calculating innundation fractions
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: real2string, int2string, logical2string
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_top_config

  TYPE, EXTENDS(t_jsb_config) :: t_top_config
    LOGICAL           :: dbg_topmodel           ! Debug topmodel processes and tasks
    REAL(wp)          :: efold                  ! topmodel efolding depth
    REAL(wp)          :: wt_min                 ! minimum water table required for inundation
    REAL(wp)          :: wt_max                 ! maximum water table for inundation (larger is considered lake/river)
    REAL(wp)          :: cti_min                ! minimum mean CTI value to consider grid cell
    REAL(wp)          :: cti_max                ! maximum CTI for inundation (larger is considered lake/river)
    REAL(wp)          :: cti_min_flooding       ! minimum CTI for floodability
    REAL(wp)          :: saturation_threshold   ! threshold above which soil layer is assumed saturated
    REAL(wp)          :: relaxation_time        ! time scale for relaxation of aet and precip [d]
    REAL(wp)          :: wetland_threshold      ! Minimum for relaxed P/E allowing wetlands []
    LOGICAL           :: init_runmean           ! initialise the running means
    LOGICAL           :: use_cti_basin          ! use basin mean CTI index in inundated area calc.
    LOGICAL           :: use_openwater          ! use open water fraction to limit inundated fraction
    LOGICAL           :: use_wt_max             ! use wt_max to limit inundated fraction
    LOGICAL           :: use_cti_max            ! use cti_max to limit inundated fraction
    LOGICAL           :: use_cti_flood          ! use cti_min_flooding to limit inundated fraction
    REAL(wp)          :: peat_runoff_exp        ! exponent in peat runoff equation
    REAL(wp)          :: peat_runoff_fac        ! factor in peat runoff equation
    LOGICAL           :: inun_prescribed        ! prescribe inundation
    LOGICAL           :: lpeat                  ! do peat calculations
    LOGICAL           :: ltime_means            ! determine time means
    REAL(wp)          :: lat_const_wetland_type ! Latitude separating types of calculations for permanent wetland
    INTEGER           :: month_summer_start     ! Starting #month for summer specific calculations
    INTEGER           :: month_summer_end       ! Ending #month for summer specific calculations
    INTEGER           :: years_long_term        ! #years to include in long term means
    REAL(wp)          :: soil_depth_min         ! Minimum soil depth to deterimine water table and inundation [m]
    REAL(wp)          :: scale_inun             ! Scale inundated fraction with this number []
  CONTAINS
    PROCEDURE :: Init => Init_topmodel_config
  END TYPE t_top_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_topconfig_class'

CONTAINS

  SUBROUTINE Init_topmodel_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input

    CLASS(t_top_config), INTENT(inout) :: config

    ! Namelist variables
    LOGICAL                     :: active
    LOGICAL                     :: dbg_topmodel
    REAL(wp)                    :: efold
    REAL(wp)                    :: wt_min
    REAL(wp)                    :: wt_max
    REAL(wp)                    :: cti_min
    REAL(wp)                    :: cti_max
    REAL(wp)                    :: cti_min_flooding
    REAL(wp)                    :: saturation_threshold
    REAL(wp)                    :: relaxation_time
    REAL(wp)                    :: wetland_threshold
    LOGICAL                     :: init_runmean
    LOGICAL                     :: use_cti_basin
    LOGICAL                     :: use_openwater
    LOGICAL                     :: use_wt_max
    LOGICAL                     :: use_cti_max
    LOGICAL                     :: use_cti_flood
    REAL(wp)                    :: peat_runoff_exp
    REAL(wp)                    :: peat_runoff_fac
    LOGICAL                     :: inun_prescribed
    LOGICAL                     :: lpeat
    LOGICAL                     :: ltime_means
    REAL(wp)                    :: lat_const_wetland_type
    INTEGER                     :: month_summer_start
    INTEGER                     :: month_summer_end
    INTEGER                     :: years_long_term
    REAL(wp)                    :: soil_depth_min
    REAL(wp)                    :: scale_inun
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_topmodel_nml/                      &
      & active, dbg_topmodel,                        &
      & efold, wt_min, wt_max,                       &
      & cti_min, cti_max, cti_min_flooding,          &
      & saturation_threshold, relaxation_time,       &
      & wetland_threshold, init_runmean,             &
      & use_cti_basin, use_openwater,                &
      & use_wt_max, use_cti_max, use_cti_flood,      &
      & peat_runoff_exp, peat_runoff_fac,            &
      & inun_prescribed, lpeat, ltime_means,         &
      & lat_const_wetland_type, years_long_term,     &
      & month_summer_start, month_summer_end,        &
      & soil_depth_min, scale_inun,                  &
      & ic_filename, bc_filename

    ! Locals
    INTEGER :: nml_handler, nml_unit, istat, i
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_topmodel_config'

    CALL message(TRIM(routine), 'Starting topmodel configuration')

    ! Set defaults
    active                 = .FALSE.
    dbg_topmodel           = .FALSE.
    cti_min                =   5.0_wp
    cti_max                =  35.0_wp
    cti_min_flooding       =   9.0_wp
    wt_min                 =   0._wp
    wt_max                 =   0.5_wp
    efold                  =   2.0_wp
    saturation_threshold   =   0.95_wp
    relaxation_time        = 365.2425_wp
    wetland_threshold      =   1.0_wp
    init_runmean           = .FALSE.
    use_cti_basin          = .FALSE.
    use_openwater          = .FALSE.
    use_wt_max             = .FALSE.
    use_cti_max            = .FALSE.
    use_cti_flood          = .TRUE.
    peat_runoff_exp        =  20.0_wp
    peat_runoff_fac        =   1.e-2_wp
    inun_prescribed        = .FALSE.
    lpeat                  = .FALSE.
    ltime_means            = .FALSE.
    lat_const_wetland_type =  40._wp
    years_long_term        =   5
    month_summer_start     =   6 ! Summer starts in June
    month_summer_end       =   9 ! Summer ends in September
    soil_depth_min         =   0.11_wp
    scale_inun             =   1._wp
    ic_filename            = 'ic_land_soil.nc'
    bc_filename            = 'bc_land_soil.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_topmodel_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_topmodel_nml)

    CALL close_nml(nml_handler)

    config%active                 = active
    config%dbg_topmodel           = dbg_topmodel
    config%cti_min                = cti_min
    config%efold                  = efold
    config%wt_min                 = wt_min
    config%wt_max                 = wt_max
    config%cti_max                = cti_max
    config%cti_min_flooding       = cti_min_flooding
    config%saturation_threshold   = saturation_threshold
    config%relaxation_time        = relaxation_time
    config%wetland_threshold      = wetland_threshold
    config%use_cti_basin          = use_cti_basin
    config%use_openwater          = use_openwater
    config%use_wt_max             = use_wt_max
    config%use_cti_max            = use_cti_max
    config%use_cti_flood          = use_cti_flood
    config%peat_runoff_exp        = peat_runoff_exp
    config%peat_runoff_fac        = peat_runoff_fac
    config%inun_prescribed        = inun_prescribed
    config%init_runmean           = init_runmean
    config%lpeat                  = lpeat
    config%ltime_means            = ltime_means
    config%lat_const_wetland_type = lat_const_wetland_type
    config%month_summer_start     = month_summer_start
    config%month_summer_end       = month_summer_end
    config%years_long_term        = years_long_term
    config%scale_inun             = scale_inun
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename

    ! Tests and messsages
    IF (active) THEN
      IF (use_cti_flood .AND. use_wt_max) THEN
        CALL message(routine,'use_flood_cti and use_wt_max are mutually exclusive -> disabling use_wt_max')
        use_wt_max = .FALSE.
      END IF
      IF (use_cti_flood .AND. use_cti_max) THEN
        CALL message(routine,'use_cti_flood and use_cti_max are mutually exclusive -> disabling use_cti_max')
        use_cti_max = .FALSE.
      END IF
      IF (use_wt_max .AND. use_cti_max) THEN
        CALL message(routine,'use_wt_max and use_cti_max are mutually exclusive -> disabling use_wt_max')
        use_wt_max = .FALSE.
      END IF
      IF (lpeat .AND. inun_prescribed) THEN
        CALL message(routine,'lpeat and prescribed inundation are mutually exclusive -> disabling lpeat')
        lpeat = .FALSE.
      END IF
      CALL message(routine,'Peatland water table submodel: lpeat '//logical2string(lpeat))
      IF (inun_prescribed) CALL message(routine,'Using prescribed inundation fraction. Many other parameters ignored')
      IF (dbg_topmodel) THEN
        IF (.NOT. inun_prescribed) THEN
          CALL message(routine,'Detailed topmodel settings: ')
          CALL message(routine,'  CTI data from         : '//ic_filename)
          CALL message(routine,'  cti_min               : '//real2string(cti_min))
          CALL message(routine,'  efold                 : '//real2string(efold))
          CALL message(routine,'  wt_min                : '//real2string(wt_min))
          CALL message(routine,'  saturation_threshold  : '//real2string(saturation_threshold))
          CALL message(routine,'  Max extent option: ')
          CALL message(routine,'    use_cti_flood   : '//logical2string(use_cti_flood))
          IF (use_cti_flood) CALL message(routine,'    cti_min_flooding: '//real2string(cti_min_flooding))
          CALL message(routine,'    use_wt_max      : '//logical2string(use_wt_max))
          IF (use_wt_max) CALL message(routine,'      wt_max        : '//real2string(wt_max))
          CALL message(routine,'    use_cti_max     : '//logical2string(use_cti_max))
          IF (use_cti_max) CALL message(routine,'      cti_max       :'//real2string(cti_max))
          CALL message(routine,'    use_openwater   : '//logical2string(use_openwater))
          CALL message(routine,'    use_cti_basin   : '//logical2string(use_cti_basin))
          IF (lpeat) THEN
            CALL message(routine,'    peat_runoff_exp   : '//real2string(peat_runoff_exp))
            CALL message(routine,'    peat_runoff_fac   : '//real2string(peat_runoff_fac))
          END IF
          CALL message(routine,'  lat_const_wetland_type: '//real2string(lat_const_wetland_type))
          CALL message(routine,'  years_long_term       : '//int2string(years_long_term))
          CALL message(routine,'  month_summer_start    : '//int2string(month_summer_start))
          CALL message(routine,'  month_summer_end      : '//int2string(month_summer_end))
          ! StW: Guess settings below are duplicated with general model settings
          CALL message(routine,'  Determine time mean diagnostics: ltime_means '//logical2string(ltime_means))
          IF (ltime_means) CALL message(routine,'  init_runmean '//logical2string(init_runmean))
        END IF
        CALL message(routine,'Topmodel configured')
      ENDIF
    ELSE
      CALL message(routine,'Topmodel is not used')
    ENDIF

    !! -----------------------------------------------------------------------------------------------------------------------------
    IF (.NOT. active) RETURN

  END SUBROUTINE Init_topmodel_config

#endif
END MODULE mo_top_config_class
