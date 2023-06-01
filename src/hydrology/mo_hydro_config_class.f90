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
  USE mo_util_string,       ONLY: tolower
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hydro_config

  TYPE, EXTENDS(t_jsb_config) :: t_hydro_config
    CHARACTER(len=filename_max) :: bc_sso_filename !< For elevation and oro_stddev
    REAL(wp)               :: &
      & w_skin_max,           & !< Maximum water holding capacity of the soil and canopy skins, as well as
                                !  maximum snow depth on leaves  [m water equivalent]
      & w_soil_limit,         & !< Upper limit for maximum soil moisture content
      & w_soil_crit_fract,    & !< Fraction of field capacity at which soil starts to dry,
                                !  i.e. at which water stress starts to be above zero
      & w_soil_wilt_fract,    & !< Fraction of field capacity at which soil is at wilting point,
                                !  i.e. when no transpiration occurs anymore
      & snow_depth_max,       & !< Limit for snow depth [m water equivalent] to avoid grid cells with
                                !  infinitely growing snow depth; -1. for no limitation
      & manning_n,            & ! Manning cefficient for surface roughness: 0.035 -> floodplains - pasture,
                                ! 0.05 -> floodplains - light brush, 0.075 -> floodplains - heavy brush, 0.15 -> floodplains - trees
      & dt_sat,               & !< Time required for hydraulic conductivity to reach max [sec]
      & max_slope,            & !< maximum slope considered for infiltration/runoff estimation
      & max_pond_depth          !< prescribed maximum surface pond depth

    CHARACTER(len=10)     ::  &
      & scheme_stom_cond        !< Scheme to use for computation of stomatal conductance (echam5 only for now)
    LOGICAL               ::  &
      & enforce_water_budget, & !< stop in case of water conservation problem
      & l_organic,            & !< Use organic component of soil layers
      & l_socmap,             & !< Replace computed organic layer with 3D map of soil carbon fraction
      & l_glac_in_soil,       & !< Treat glacier as part of soil, not separate glacier tile
      & l_infil_subzero,      & !< Allow infiltration at soil temperatures below 0degC
      & l_ponds,              & !< Additional surface water storage
      & l_jsb_wet_test          !< Create idealized, wet jsb-setup

    INTEGER               ::  &
      & runoff_scheme,        & !< Choice of runoff and infiltration scheme: 'ARNO' (default), 'QArctic' (only with hydrotiles)
      & soilhydmodel,         & !< Choice of soil hydrological model: 'BrooksCorey' (default), 'ClappHornberger' or 'VanGenuchten'
      & interpol_mean           !< Interpolation scheme, to optain K & D at layer interface: 'upstream' (default) or 'arithmetic' mean

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
    USE mo_hydro_constants,    ONLY: BrooksCorey_, ClappHornberger_, VanGenuchten_, &
      &                              Upstream_, Arithmetic_, Arno_, QArctic_

    CLASS(t_hydro_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename, bc_sso_filename
    REAL(wp)                    :: w_skin_max, snow_depth_max
    REAL(wp)                    :: manning_n, dt_sat, max_slope, max_pond_depth
    REAL(wp)                    :: w_soil_limit, w_soil_crit_fract, w_soil_wilt_fract
    CHARACTER(len=15)           :: scheme_stom_cond, runoff_scheme, soilhydmodel, interpol_mean
    LOGICAL                     :: enforce_water_budget, l_organic, l_socmap, l_glac_in_soil, l_infil_subzero, &
      &                            l_ponds, l_jsb_wet_test

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
      & manning_n,                  &
      & dt_sat,                     &
      & max_slope,                  &
      & max_pond_depth,             &
      & scheme_stom_cond,           &
      & enforce_water_budget,       &
      & l_organic,                  &
      & l_socmap,                   &
      & l_glac_in_soil,             &
      & l_infil_subzero,            &
      & l_ponds,                    &
      & l_jsb_wet_test,             &
      & runoff_scheme,              &
      & soilhydmodel,               &
      & interpol_mean

    INTEGER :: nml_handler, nml_unit, istat
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)
    REAL(wp), ALLOCATABLE :: depths(:), mids(:)
    INTEGER :: nsoil

    TYPE(t_jsb_vgrid), POINTER :: soil_w

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hydro_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting hydro configuration')

    ! Set defaults
    active                 = .TRUE.
    bc_filename            = 'bc_land_hydro.nc'  !R: gibts im Moment noch nicht
    ic_filename            = 'ic_land_hydro.nc'  !R: gibts im Moment noch nicht
    bc_sso_filename        = 'bc_land_sso.nc'
    w_skin_max             = 2.E-4_wp            ! 0.2 mm corresponds to Roesch et al. 2001
    w_soil_limit           = -1.0_wp
    w_soil_crit_fract      = 0.75_wp
    w_soil_wilt_fract      = 0.35_wp
    snow_depth_max         = -1.0_wp             ! -1: no limitation
    manning_n              = 0.075_wp
    dt_sat                 = 600._wp             ! [sec]
    max_slope              = 0.45_wp             ! [m/m]
    max_pond_depth         = 0.1_wp              ! [m]
    scheme_stom_cond       = 'echam5'
    enforce_water_budget   = config%model_config%enforce_water_budget
    l_organic              = .TRUE.
    l_socmap               = .TRUE.
    l_glac_in_soil         = .FALSE.
    l_infil_subzero        = .TRUE.
    l_ponds                = .FALSE.
    l_jsb_wet_test         = .FALSE.
    runoff_scheme          = "Arno"
    soilhydmodel           = "BrooksCorey"
    interpol_mean          = "upstream"

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_hydro_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_hydro_nml)

    IF (TRIM(config%model_config%usecase) == 'jsbach_qarctic') THEN
      ! Q-Arctic use case needs different defaults for some hydro parameters
      l_ponds              = .TRUE.
      runoff_scheme        = 'QArctic'
    ENDIF

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%bc_sso_filename     = bc_sso_filename

    config%w_skin_max          = w_skin_max
    CALL message(TRIM(routine), 'Maximum water holding capacity of soil or canopy skin: '//TRIM(real2string(w_skin_max)))

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

    config%manning_n      = manning_n
    CALL message(TRIM(routine), 'Manning coefficient for lateral flow: '//TRIM(real2string(manning_n)))

    config%dt_sat         = dt_sat
    CALL message(TRIM(routine), 'Time required for hydraulic conductivity to reach max [sec]: '//TRIM(real2string(dt_sat)))

    config%max_slope      = max_slope
    CALL message(TRIM(routine), 'Maximum slope considered in runoff and infiltration computation [sec]: '&
      & //TRIM(real2string(max_slope)))

    config%max_pond_depth = max_pond_depth
    CALL message(TRIM(routine), 'Prescribed maximum pond depth [m]: '//TRIM(real2string(max_pond_depth)))

    config%scheme_stom_cond       = scheme_stom_cond

    config%enforce_water_budget   = enforce_water_budget
    IF (.NOT. enforce_water_budget) THEN
      CALL message(TRIM(routine), 'WARNING: Model will not stop due to land surface water balance violation.')
    END IF

    config%l_organic              = l_organic
    IF (l_organic) THEN
      CALL message(TRIM(routine), 'Using organic component in soil hydrology')
    ELSE
      IF (l_socmap) CALL message(TRIM(routine), 'Setting l_socmap=false since l_organic=false in soil hydrology')
      l_socmap = .FALSE.
    END IF

    config%l_socmap = l_socmap
    IF (l_socmap .AND. l_organic) &
      & CALL message(TRIM(routine), 'Replacing organic carbon fractions with map from bc file: ' //TRIM(logical2string(l_socmap)))

    config%l_glac_in_soil       = l_glac_in_soil
    config%l_infil_subzero      = l_infil_subzero
    IF (l_infil_subzero) CALL message(TRIM(routine), 'Allow infiltration at sub-zero temperatures')

    config%l_ponds = l_ponds

    IF (l_ponds) THEN
      CALL message(TRIM(routine), '*** WARNING: Use of ponds is still experimental and not tested!!')
      config%max_pond_depth = max_pond_depth
      CALL message(TRIM(routine), 'Prescribed maximum pond depth [m]: '//TRIM(real2string(config%max_pond_depth)))
    ELSE
      config%max_pond_depth = 0._wp
      CALL message(TRIM(routine), 'Experimental pond scheme is disabled.')
    END IF

    config%l_jsb_wet_test = l_jsb_wet_test
    IF (l_jsb_wet_test) THEN
      CALL message(TRIM(routine), '*** This is a WET run***')
    END IF

    SELECT CASE(TRIM(tolower(runoff_scheme)))
    CASE("arno")
      config%runoff_scheme = Arno_
      CALL message(TRIM(routine), 'Soil hydrology uses ARNO scheme for infiltration and surface runoff compuation')
    CASE("qarctic")
      config%runoff_scheme = QArctic_
      CALL message(TRIM(routine), 'Soil hydrology uses QArctic scheme for infiltration and surface runoff compuation')
    CASE default
      CALL finish(TRIM(routine), 'Runoff scheme '//TRIM(runoff_scheme)//' not implemented')
    END SELECT

    SELECT CASE(TRIM(tolower(soilhydmodel)))
    CASE("brookscorey")
      config%soilhydmodel = BrooksCorey_
      CALL message(TRIM(routine), 'Soil hydrology uses the Brooks & Corey model')
    CASE("clapphornberger")
      config%soilhydmodel = ClappHornberger_
      CALL message(TRIM(routine), 'Soil hydrology uses the Clapp & Hornberger model')
    CASE("vangenuchten")
      config%soilhydmodel = VanGenuchten_
      CALL message(TRIM(routine), 'Soil hydrology uses the Van Genuchten model')
    CASE default
      CALL finish(TRIM(routine), 'Soil hydrology model '//TRIM(soilhydmodel)//' not implemented')
    END SELECT

    SELECT CASE(TRIM(tolower(interpol_mean)))
    CASE("upstream")
      config%interpol_mean = Upstream_
      CALL message(TRIM(routine), 'Soil diffusivity/conductivity interpolation uses upstream means')
    CASE("arithmetic")
      config%interpol_mean = Arithmetic_
      CALL message(TRIM(routine), 'Soil diffusivity/conductivity interpolation uses arithmetic means')
    CASE default
      CALL finish(TRIM(routine), 'Soil diffusivity/conductivity interpolation method '//TRIM(interpol_mean)//' not implemented')
    END SELECT

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
