!>
!! @brief Contains structures and methods for hydrology config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-30)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hydro_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_util,              ONLY: real2string, logical2string
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hydro_config

  TYPE, EXTENDS(t_jsb_config) :: t_hydro_config
    CHARACTER(len=filename_max) :: bc_sso_filename !< For elevation and oro_stddev
    REAL(wp)            :: &
      & w_skin_max,        & !< Water holding capacity of ground skin reservoir [m water equivalent]
      & w_soil_limit,      & !< Upper limit for maximum soil moisture content
      & w_soil_crit_fract, & !< Fraction of field capacity at which soil starts to dry,
                             !  i.e. at which water stress starts to be above zero
      & w_soil_wilt_fract, & !< Fraction of field capacity at which soil is at wilting point,
                             !  i.e. when no transpiration occurs anymore
      & snow_depth_max       !< Limit for snow depth [m water equivalent] to avoid grid cells with
                             !  infinitely growing snow depth; -1. for no limitation
    CHARACTER(len=10)     :: &
      & scheme_stom_cond     !< Scheme to use for computation of stomatal conductance (echam5 only for now)
    LOGICAL               :: &
      & l_organic,           & !< Use organic component of soil layers
      & l_glac_in_soil,      & !< Treat glacier as part of soil, not separate glacier tile
!pdv >
      & l_hyd_compat420,        & !< Reference version: JSBACH 4.20, i.e. without KoPF changes
      & l_runoff_with_ice,      & !< Include soil ice in runoff and infiltration calculations
      & l_ws_max_depend_on_org, & !< ws_max re-calculation (in calc_soil_hydrology) depends on organic soil fraction
      & l_correct_ws_from_top,  & !< Start soil water correction from top layer (instead of equal reduction from all layers)
      & l_ice_impedance,        & !< Account for ice impedance
      & l_diffuse_all_ws,       & !< Calculate soil water diffusion from all soil water, including immobile
      & l_infil_subzero,        & !< Allow infiltration at soil temperatures below 0°C
      & l_discount_scool,       & !< supercooled water is neither available to plants nor percolates
      & l_drain_lat,            & !< Lateral drainage, i.e. drainage from upper to lower soil layers;
                                  !  if false, only the lowest soil level drains.
      & l_stress_from_all_ws,   & !< Calculate water stress from all soil water, including water below
                                  !  the root zone and including supercooled water.
      & l_percolation_theory      !< The saturated hydraulic conductivity is calculated following
                                  !  the percolation theory
!pdv <
  CONTAINS
    PROCEDURE             :: Init => Init_hydro_config
  END type t_hydro_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_config_class'

CONTAINS

  SUBROUTINE Init_hydro_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input

    CLASS(t_hydro_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename, bc_sso_filename
    REAL(wp)                    :: w_skin_max, snow_depth_max
    REAL(wp)                    :: w_soil_limit, w_soil_crit_fract, w_soil_wilt_fract
    CHARACTER(len=10)           :: scheme_stom_cond
    LOGICAL                     :: l_organic, l_glac_in_soil
!pdv >
    LOGICAL     :: l_hyd_compat420
    LOGICAL     :: l_infil_subzero, l_runoff_with_ice, l_ws_max_depend_on_org, l_correct_ws_from_top, l_ice_impedance
    LOGICAL     :: l_diffuse_all_ws, l_discount_scool, l_drain_lat, l_stress_from_all_ws, l_percolation_theory
!pdv <

    NAMELIST /jsb_hydro_nml/        &
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & bc_sso_filename,            &
      & w_skin_max,                 &
      & w_soil_limit,               &
      & w_soil_crit_fract,          &
      & w_soil_wilt_fract,          &
      & snow_depth_max,             &
      & scheme_stom_cond,           &
      & l_organic,                  &
      & l_glac_in_soil,             &
!pdv >
      & l_hyd_compat420,            &
      & l_infil_subzero,            &
      & l_runoff_with_ice,          &
      & l_ws_max_depend_on_org,     &
      & l_correct_ws_from_top,      &
      & l_ice_impedance,            &
      & l_diffuse_all_ws,           &
      & l_discount_scool,           &
      & l_drain_lat,                &
      & l_stress_from_all_ws,       &
      & l_percolation_theory
!pdv <

    INTEGER :: nml_handler, nml_unit, istat
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)
    REAL(wp), ALLOCATABLE :: depths(:), mids(:)
    INTEGER :: nsoil

    TYPE(t_jsb_vgrid), POINTER :: soil_w

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hydro_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting hydro configuration')

    ! Set defaults
    active            = .TRUE.
    bc_filename       = 'bc_land_hydro.nc'  !R: gibts im Moment noch nicht
    ic_filename       = 'ic_land_hydro.nc'  !R: gibts im Moment noch nicht
    bc_sso_filename   = 'bc_land_sso.nc'
    w_skin_max        = 2.E-4_wp
    w_soil_limit      = -1.0_wp
    w_soil_crit_fract = 0.75_wp
    w_soil_wilt_fract = 0.35_wp
    snow_depth_max    = -1.0_wp             ! -1: no limitation
    scheme_stom_cond  = 'echam5'
    l_organic         = .FALSE.
    l_glac_in_soil    = .FALSE.
!pdv >
    l_hyd_compat420        = .TRUE.
    l_infil_subzero        = .FALSE.
    l_runoff_with_ice      = .TRUE.
    l_ws_max_depend_on_org = .TRUE.
    l_correct_ws_from_top  = .FALSE.
    l_ice_impedance        = .FALSE.
    l_diffuse_all_ws       = .TRUE.
    l_discount_scool       = .FALSE.
    l_drain_lat            = .TRUE.
    l_stress_from_all_ws   = .TRUE.
    l_percolation_theory   = .FALSE.
!pdv <

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_hydro_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_hydro_nml)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%bc_sso_filename     = bc_sso_filename

    config%w_skin_max          = w_skin_max
    CALL message(TRIM(routine), 'Maximum content of soil skin reservoir: '//TRIM(real2string(w_skin_max)))

    config%w_soil_limit        = w_soil_limit
    CALL message(TRIM(routine), 'Upper limit for maximum soil moisture content: '//TRIM(real2string(w_soil_limit)))

    config%w_soil_crit_fract   = w_soil_crit_fract
    CALL message(TRIM(routine), 'Fraction of field capacity at critical point: '//TRIM(real2string(w_soil_crit_fract)))

    config%w_soil_wilt_fract   = w_soil_wilt_fract
    CALL message(TRIM(routine), 'Fraction of field capacity at permanent wilting point: '//TRIM(real2string(w_soil_wilt_fract)))

    config%snow_depth_max      = snow_depth_max
    IF (snow_depth_max >= 0._wp) THEN
      CALL message(TRIM(routine), 'Snow depth limitation: '//TRIM(real2string(snow_depth_max)))
    END IF

    config%scheme_stom_cond    = scheme_stom_cond

    config%l_organic           = l_organic
    IF (l_organic) CALL message(TRIM(routine), 'Using organic component in soil hydrology')

    config%l_glac_in_soil      = l_glac_in_soil
!pdv >
    config%l_hyd_compat420     = l_hyd_compat420
    IF (l_hyd_compat420) THEN
      CALL message(TRIM(routine), 'Hydrology settings as in JSBACH 4.20')
      IF (l_infil_subzero) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_infil_subzero=F. Namelist definition of l_infil_subzero ignored.')
        l_infil_subzero        = .FALSE.
      END IF
      IF (.NOT. l_runoff_with_ice) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_runoff_with_ice=T. Namelist definition of l_runoff_with_ice ignored.')
        l_runoff_with_ice      = .TRUE.
      END IF
      IF (.NOT. l_ws_max_depend_on_org) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_ws_max_depend_on_org=T. Namelist definition of l_ws_max_depend_on_org ignored.')
        l_ws_max_depend_on_org = .TRUE.
      END IF
      IF (l_correct_ws_from_top) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_correct_ws_from_top=F. Namelist definition of l_correct_ws_from_top ignored.')
        l_correct_ws_from_top  = .FALSE.
      END IF
      IF (l_ice_impedance) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_ice_impedance=F. Namelist definition of l_ice_impedance ignored.')
        l_ice_impedance        = .FALSE.
      END IF
      IF (.NOT. l_diffuse_all_ws) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_diffuse_all_ws=T. Namelist definition of l_diffuse_all_ws ignored.')
        l_diffuse_all_ws       = .TRUE.
      END IF
      IF (l_discount_scool) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_discount_scool=F. Namelist definition of l_discount_scool ignored.')
        l_discount_scool       = .FALSE.
      END IF
      IF (.NOT. l_drain_lat) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_drain_lat=T. Namelist definition of l_drain_lat ignored.')
        l_drain_lat            = .TRUE.
      END IF
      IF (l_percolation_theory) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_percolation_theory=F. Namelist definition of l_percolation_theory ignored.')
        l_percolation_theory   = .FALSE.
      END IF
      IF (.NOT. l_stress_from_all_ws) THEN
        CALL message(TRIM(routine), &
         & 'WARNING: with l_hyd_compat420  l_stress_from_all_ws=T. Namelist definition of l_stress_from_all_ws ignored.')
        l_stress_from_all_ws   = .TRUE.
      END IF
    END IF
    config%l_infil_subzero        = l_infil_subzero
    CALL message(TRIM(routine), &
     & 'Infiltration possible at soil temperatures below 0 degC: '//TRIM(logical2string(l_infil_subzero)))
    config%l_runoff_with_ice      = l_runoff_with_ice
    CALL message(TRIM(routine), &
     & 'Calculate runoff including ice: '//TRIM(logical2string(l_runoff_with_ice)))
    config%l_ws_max_depend_on_org = l_ws_max_depend_on_org
    CALL message(TRIM(routine), &
     & 'Re-calculate ws_max from initialfile dep. on organic: '//TRIM(logical2string(l_ws_max_depend_on_org)))
    config%l_correct_ws_from_top  = l_correct_ws_from_top
    CALL message(TRIM(routine), &
     & 'Start correction of soil moisture from top level: '//TRIM(logical2string(l_correct_ws_from_top)))
    config%l_ice_impedance        = l_ice_impedance
    CALL message(TRIM(routine), &
     & 'Account for ice impedance: '//TRIM(logical2string(l_ice_impedance)))
    config%l_diffuse_all_ws       = l_diffuse_all_ws
    CALL message(TRIM(routine), &
     & 'Calculate diffusion of all water, including supercooled water: '//TRIM(logical2string(l_diffuse_all_ws)))
    config%l_discount_scool       = l_discount_scool
    CALL message(TRIM(routine), &
     & 'Account for supercooled water: '//TRIM(logical2string(l_discount_scool)))
    config%l_drain_lat            = l_drain_lat
    CALL message(TRIM(routine), &
     & 'Lateral drainage is possible: '//TRIM(logical2string(l_drain_lat)))
    config%l_percolation_theory   = l_percolation_theory
    CALL message(TRIM(routine), &
     & 'Calculate percolation from Percolation Theory: '//TRIM(logical2string(l_percolation_theory)))
    config%l_stress_from_all_ws   = l_stress_from_all_ws
    CALL message(TRIM(routine), &
     & 'Account for all water for water stress calculation, including water below the rootzone and supercooled water: ' &
     & //TRIM(logical2string(l_stress_from_all_ws)))
!pdv <

    CALL close_nml(nml_handler)

    IF (.NOT. active) RETURN

    input_file = jsb_netcdf_open_input(ic_filename)

    ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
    nsoil = SIZE(ptr_1D)

    IF (l_organic .AND. nsoil < 3) &
      & CALL finish(TRIM(routine), 'At least three soil layers required for using organic soil component')

    ALLOCATE(depths(nsoil+1))
    ALLOCATE(mids(nsoil))
    depths(1) = 0._wp
    depths(2:nsoil+1) = ptr_1D(1:nsoil)
    mids(1:nsoil) = (depths(1:nsoil) + depths(2:nsoil+1)) / 2._wp
    DEALLOCATE(ptr_1D)

    CALL input_file%Close()

    soil_w  => new_vgrid('soil_depth_water', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
      & levels  = mids                 (1:nsoil  ),                         &
      & lbounds = depths               (1:nsoil  ),                         &
      & ubounds = depths               (2:nsoil+1),                         &
      & units='m')
    CALL register_vgrid(soil_w)
    WRITE(message_text, *) 'Soil levels in hydrology (upper) [m]: ', soil_w%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in hydrology (mid)   [m]: ', soil_w%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in hydrology (lower) [m]: ', soil_w%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil level depths in hydrology   [m]: ', soil_w%dz
    CALL message(TRIM(routine), message_text)
    DEALLOCATE(depths, mids)

  END SUBROUTINE Init_hydro_config

#endif
END MODULE mo_hydro_config_class
