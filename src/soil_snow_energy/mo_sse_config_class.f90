!>
!! @brief Contains structures and methods for soil and snow energy config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_sse_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: real2string, int2string
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sse_config

  TYPE, EXTENDS(t_jsb_config) :: t_sse_config
    INTEGER  :: nsnow           !< Maximum number of snow layers
    LOGICAL  :: l_snow          !< T: Consider snow layers in soil thermal calculations
    ! l_dynsnow, l_heat_cap_dyn and l_heat_cond_dyn are only effective for l_snow=.true.
    LOGICAL  :: l_dynsnow       !< T: Calculate snow parameters heat cond. and heat cap. dynamically (depending on snow density)
    LOGICAL  :: l_heat_cap_dyn  !< T: Use dynamic calculation of soil heat capacity, F: use mineral soil map
    LOGICAL  :: l_heat_cond_dyn !< T: Use dynamic calculation of soil heat conductivity, F: use mineral soil map
    !
    LOGICAL  :: l_heat_cap_map  !< T: Use soil heat capacity from input map, F: derive from FAO
    LOGICAL  :: l_heat_cond_map !< T: USe soil heat conductivity from input map, F: derive from heat cap. and FAO thermal diff.
    LOGICAL  :: l_freeze        !< T: Consider freezing and thawing in thermal soil calculations
    LOGICAL  :: l_supercool     !< T: Allow for supercooled soil water
    REAL(wp) :: w_soil_critical !< Critical water/ice content in upper soil layer for correction of new surface
                                !! temperature for freezing/melting [m water equivalent]
    !! quincy
    ! INTEGER       :: nsoil_energy         !< number of soil layers for soil energy calculations   ! depreciated feature, but needed for now, see quincy-trac ticket #521
    ! INTEGER       :: nsoil_water          !< number of soil layers for soil moisture calculations ! depreciated feature, but needed for now, see quincy-trac ticket #521
    REAL(wp)      :: soil_depth           !< actual soil and rooting depth
    REAL(wp)      :: wsr_capacity         !< Water holding capacity of ground skin reservoir [m water equivalent]
    REAL(wp)      :: wsn_capacity         !< Water holding capacity of ground snow reservoir (max snow depth) [m water equivalent]
    REAL(wp)      :: soil_awc_prescribe, &
                     soil_theta_prescribe
    REAL(wp)      :: soil_sand            !< soil sand proportion - site specific from forcing info file
    REAL(wp)      :: soil_silt            !< soil silt proportion - site specific from forcing info file
    REAL(wp)      :: soil_clay            !< soil clay proportion - recalculated from: clay = 1.0 -sand -silt
    REAL(wp)      :: bulk_density    

  CONTAINS
    PROCEDURE :: Init => Init_sse_config
  END TYPE t_sse_config

  INTEGER, PARAMETER :: max_snow_layers = 20
  INTEGER, PARAMETER :: max_soil_layers = 20

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_config_class'

CONTAINS

  SUBROUTINE Init_sse_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND, ZAXIS_GENERIC
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_sse_constants,      ONLY: snow_depth_min
    USE mo_jsb_math_constants, ONLY: eps8

    CLASS(t_sse_config), INTENT(inout) :: config

    LOGICAL  :: active, l_snow, l_dynsnow, l_heat_cap_dyn, l_heat_cond_dyn, l_heat_cap_map, l_heat_cond_map, l_freeze, l_supercool
    REAL(wp) :: w_soil_critical
    INTEGER  :: nsnow
    REAL(wp) :: dz_snow(max_snow_layers)
    !! quincy
    ! INTEGER                     :: nsoil_energy    , &
    !                                nsoil_water
    ! REAL(wp)                    :: dz_energy(max_soil_layers), &                        ! vector of soil layer thicknesses 
    !                                dz_water(max_soil_layers)
    ! REAL(wp)                    :: soil_layer_profile_ubound_estimate(max_soil_layers)  ! intial estimate of soil profile using the upper bound of each soil layer [m]
    ! REAL(wp), ALLOCATABLE       :: ubounds_soil_lay_energy(:), &                        ! layer upper bound (note! upper bound > lower bound)
    !                                ubounds_soil_lay_water(:)
    REAL(wp)                    :: wsr_capacity    , &
                                   wsn_capacity    , &
                                   soil_depth   
    REAL(wp)                    :: soil_awc_prescribe, &
                                   soil_theta_prescribe
    REAL(wp)                    :: soil_sand       , &
                                   soil_silt       , &
                                   soil_clay
    REAL(wp)                    :: bulk_density
    INTEGER                     :: isoil   ! looping
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    REAL(wp)                    :: k_soil_profile       ! for soil layer calculations
    REAL(wp)                    :: min_layer_depth      ! for soil layer calculations
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_sse_nml/                      &
      & active,                                 &
      & nsnow, dz_snow, l_snow, l_dynsnow,      &
      & l_heat_cap_dyn, l_heat_cond_dyn,        &
      & l_heat_cap_map, l_heat_cond_map,        &
      & l_freeze, l_supercool, w_soil_critical, &
      & ic_filename, bc_filename

    !! quincy
    ! NAMELIST /sse_ctl/       & 
    !                       nsoil_energy            , &
    !                       nsoil_water             , &
    !                       dz_energy               , &
    !                       dz_water                , &
    !                       wsr_capacity            , &
    !                       wsn_capacity            , &
    !                       soil_depth              , &
    !                       soil_awc_prescribe      , &
    !                       soil_theta_prescribe    , &
    !                       soil_sand               , &
    !                       soil_silt               , &
    !                       soil_clay               , &
    !                       bulk_density       

    INTEGER :: nml_handler, nml_unit, istat, i
    TYPE(t_input_file) :: input_file
    REAL(wp), POINTER :: ptr_1D(:)
    REAL(wp), ALLOCATABLE :: depths(:), mids(:)
    INTEGER :: nsoil

    TYPE(t_jsb_vgrid), POINTER :: soil_vgrid, snow_vgrid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb             ! quincy (using values from jsb4 soil_vgrid)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_sse_config'

    ! quincy variabales for soil layer calculations
    k_soil_profile  = 0.25_wp
    min_layer_depth = 0.065_wp

    IF (debug_on()) CALL message(TRIM(routine), 'Starting soil and snow energy configuration')

    ! Set defaults
    active          = .TRUE.
    nsnow            = 5
    dz_snow(:)       = 0._wp
    dz_snow(1:nsnow) = 0.05_wp
    l_snow           = .TRUE.
    l_dynsnow        = .TRUE.
    l_heat_cap_dyn   = .TRUE.
    l_heat_cond_dyn  = .TRUE.
    l_heat_cap_map   = .FALSE.
    l_heat_cond_map  = .FALSE.
    l_freeze         = .FALSE.
    l_supercool      = .FALSE.
    w_soil_critical  = 5.85036E-3_wp
    ic_filename      = 'ic_land_soil.nc'
    bc_filename      = 'bc_land_soil.nc'
    !! quincy
    ! nsoil_energy                  = 5
    ! nsoil_water                   = 5
    ! dz_energy(1:5)                = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)
    ! dz_energy(6:max_soil_layers)  = 0._wp
    ! dz_water (1:5)                = (/0.065_wp,0.254_wp,0.913_wp,2.902_wp,5.700_wp/)
    ! dz_water (6:max_soil_layers)  = 0._wp 
    wsr_capacity                  = 2.E-4_wp
    wsn_capacity                  = 2.E-4_wp
    soil_depth                    = 5.7_wp 
    soil_awc_prescribe            = 300.0_wp
    soil_theta_prescribe          = 1.0_wp 
    soil_sand                     = 0.3_wp
    soil_silt                     = 0.4_wp
    soil_clay                     = 0.3_wp
    bulk_density                  = 1500.0_wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_sse_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_sse_nml)

    CALL close_nml(nml_handler)

    config%active        = active
    config%ic_filename   = ic_filename
    config%bc_filename   = bc_filename

    config%nsnow         = nsnow
    IF (l_snow .AND. nsnow > max_snow_layers) THEN
      CALL finish(TRIM(routine), 'Too many snow layers, maximum is '//TRIM(int2string(max_snow_layers)))
    END IF
    IF (l_snow .AND. nsnow < 3) THEN
      CALL finish(TRIM(routine), 'Number of snow layers must be larger than 2.')
    END IF
    config%l_snow        = l_snow
    IF (l_snow) THEN
      CALL message(TRIM(routine), 'Using '//TRIM(int2string(nsnow))//' snow layers.')
    ELSE
      CALL message(TRIM(routine), 'Not using multi-layer snow model')
    END IF

    IF (.NOT. l_snow .AND. l_dynsnow) THEN
      config%l_dynsnow = .FALSE.
      CALL message(TRIM(routine), 'l_dynsnow not in effect for l_snow=.false.')
    ELSE
      config%l_dynsnow       = l_dynsnow
    END IF
    IF (.NOT. l_snow .AND. l_heat_cap_dyn) THEN
      config%l_heat_cap_dyn = .FALSE.
      CALL message(TRIM(routine), 'l_heat_cap_dyn not in effect for l_snow=.false.')
    ELSE
      config%l_heat_cap_dyn  = l_heat_cap_dyn
    END IF
    IF (.NOT. l_snow .AND. l_heat_cond_dyn) THEN
      config%l_heat_cond_dyn = .FALSE.
      CALL message(TRIM(routine), 'l_heat_cond_dyn not in effect for l_snow=.false.')
    ELSE
      config%l_heat_cond_dyn = l_heat_cond_dyn
    END IF
    config%l_heat_cap_map  = l_heat_cap_map
    config%l_heat_cond_map = l_heat_cond_map
    config%l_freeze        = l_freeze
    config%l_supercool     = l_supercool

    config%w_soil_critical     = w_soil_critical
    CALL message(TRIM(routine), 'Critical water/ice content in upper soil layer for correction of '// &
      &                         'surface temperature for freezing/melting: '//TRIM(real2string(w_soil_critical)))

    !! quincy
    ! config%nsoil_energy              = nsoil_energy     ! depreciated feature, but needed for now, see ticket #521
    ! config%nsoil_water               = nsoil_water      ! depreciated feature, but needed for now, see ticket #521
    config%soil_awc_prescribe        = soil_awc_prescribe 
    config%soil_theta_prescribe      = soil_theta_prescribe 
    config%wsr_capacity              = wsr_capacity
    config%wsn_capacity              = wsn_capacity
    config%soil_depth                = soil_depth 
    config%soil_sand                 = soil_sand
    config%soil_silt                 = soil_silt
    !config%soil_clay                = not using the input value, as soil_clay must be calculated from soil_sand & soil_silt to guarantee: 1=sand+silt+clay
    config%soil_clay                 = 1.0_wp - soil_sand - soil_silt
    config%bulk_density              = bulk_density

    ! check for an error in clay content of soil
    IF(config%soil_clay < eps8) THEN
       WRITE(message_text,'(a)') 'invalid proportion of soil_clay (value < eps8)'
       CALL finish(routine, message_text)
    END IF
    ! ! check initialisation of number of soil layers
    ! IF (nsoil_water .LT. 5 .OR. nsoil_water .GT. 20) & 
    !     CALL finish(routine,'Number of soil layers needs to be a number from the range 5 to 20.')
    ! IF (nsoil_water .NE. nsoil_energy) & 
    !     CALL finish(routine,'Init with different numbers of soil layers for energy and water not yet supported for QUINCY.')


    !! -----------------------------------------------------------------------------------------------------------------------------
    IF (.NOT. active) RETURN

    input_file = jsb_netcdf_open_input(ic_filename)

    ! @todo: At the moment, the soil layers for the energy calculations are the same as for the hydrology
    ptr_1D => input_file%Read_1d(variable_name='soillev')    ! Depth of layer bottom
    nsoil = SIZE(ptr_1D)
    ALLOCATE(depths(nsoil+1))
    ALLOCATE(mids(nsoil))
    depths(1) = 0._wp
    depths(2:nsoil+1) = ptr_1D(1:nsoil)
    mids(1:nsoil) = (depths(1:nsoil) + depths(2:nsoil+1)) / 2._wp
    DEALLOCATE(ptr_1D)

    CALL input_file%Close()

    soil_vgrid  => new_vgrid('soil_depth_energy', ZAXIS_DEPTH_BELOW_LAND, nsoil, &
      & levels  = mids                 (1:nsoil  ),                              &
      & lbounds = depths               (1:nsoil  ),                              &
      & ubounds = depths               (2:nsoil+1),                              &
      & units='m')
    CALL register_vgrid(soil_vgrid)
    WRITE(message_text, *) 'Soil levels in soil energy (upper) [m]: ', soil_vgrid%lbounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in soil energy (mid)   [m]: ', soil_vgrid%levels
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil levels in soil energy (lower) [m]: ', soil_vgrid%ubounds
    CALL message(TRIM(routine), message_text)
    WRITE(message_text, *) 'Soil level depths in soil energy   [m]: ', soil_vgrid%dz
    CALL message(TRIM(routine), message_text)
    DEALLOCATE(depths, mids)

    IF (l_snow) THEN
      IF (dz_snow(1) < snow_depth_min)                                                                 &
        & CALL finish(TRIM(routine), 'Depth of first snow layer should be larger than snow_depth_min=' &
        &                            //real2string(snow_depth_min))
      ALLOCATE(depths(nsnow+1))
      ALLOCATE(mids(nsnow))
      depths(1) = 0._wp
      DO i=1,nsnow
        depths(i+1) = depths(i) + dz_snow(i)
      END DO
      mids(1:nsnow) = (depths(1:nsnow) + depths(2:nsnow+1)) / 2._wp

      snow_vgrid  => new_vgrid('snow_depth_energy', ZAXIS_GENERIC, nsnow,     &
        & levels  = mids                 (1:nsnow  ),                         &
        & lbounds = depths               (1:nsnow  ),                         &
        & ubounds = depths               (2:nsnow+1),                         &
        & units='m')
      CALL register_vgrid(snow_vgrid)
      WRITE(message_text, *) 'Snow levels in soil energy (upper) [m]: ', snow_vgrid%lbounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow levels in soil energy (mid)   [m]: ', snow_vgrid%levels
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow levels in soil energy (lower) [m]: ', snow_vgrid%ubounds
      CALL message(TRIM(routine), message_text)
      WRITE(message_text, *) 'Snow level depths in soil energy   [m]: ', snow_vgrid%dz
      CALL message(TRIM(routine), message_text)
      DEALLOCATE(depths, mids)
    END IF

    !< quincy vgrid
    !! these are "infrastructure" values, ubounds and dz must be corrected/limited by site-specific bedrock depth "soil_depth" !

    ! !! quincy vgrid with != 5 layers 
    ! ! allocate upper-bound soil-layer vectors with dimensions defined in the namelist
    ! ALLOCATE(ubounds_soil_lay_energy(nsoil_energy))
    ! ALLOCATE(ubounds_soil_lay_water(nsoil_water))
    ! 
    ! ! If not using 5 layers, estimate soil profile (i.e., dz == layer thicknesses) to mimick standard JSB4 profile
    ! !  assuming dz_water(:) = dz_energy(:)
    ! IF(nsoil_water > 5) THEN
    !   ! initial calculation of soil profile (for soil moisture calculations) calculating the upper bound of each layer
    !   DO isoil=1,nsoil_water
    !     soil_layer_profile_ubound_estimate(isoil) = min_layer_depth * & 
    !                                                 exp(k_soil_profile * isoil * REAL(max_soil_layers)/REAL(nsoil_water)) - & 
    !                                                 min_layer_depth
    !   ENDDO 
    !   ! final calculation of soil profile, calculating the layer thickness (dz(:)) and correcting for minimum soil-layer thickness
    !   dz_water(1) = min_layer_depth
    !   DO isoil=2,nsoil_water
    !     dz_water(isoil) = MAX(min_layer_depth, &
    !                            soil_layer_profile_ubound_estimate(isoil) - soil_layer_profile_ubound_estimate(isoil-1))
    !   ENDDO
    !   ! apply the same soil profile to soil profile for soil energy calculations
    !   dz_energy(:) = dz_water(:)
    ! ENDIF
    ! 
    ! ! run check for soil depth namelist-definition; may not be less/equal zero or larger than covered by soil profile
    ! IF (soil_depth .LE. 0.0_wp) CALL finish('Init_sse_config','Initialization with non positive soil depth')
    ! IF (soil_depth .GT. SUM(dz_energy)) CALL finish('Init_sse_config','Initialization of soil depth larger than soil profile not possible')
    ! 
    ! ! calc upper bound of soil layers (water and energy) from layer thickness
    ! ! note: upper bound is the larger value compared to lower bound of the same layer!
    ! ubounds_soil_lay_water(1)        = dz_water(1)
    ! DO isoil=2,nsoil_water
    !   ubounds_soil_lay_water(isoil)  = ubounds_soil_lay_water(isoil-1) + dz_water(isoil)
    ! END DO
    ! ubounds_soil_lay_energy(1)       = dz_energy(1)
    ! DO isoil=2,nsoil_energy
    !   ubounds_soil_lay_energy(isoil) = ubounds_soil_lay_energy(isoil-1) + dz_energy(isoil)
    ! END DO


    !! use values from JSBACH4 soil_vgrid init above (5 soil layers)
    !! 
    !vgrid_soil_sb => new_vgrid('soil_layer_sb', ZAXIS_GENERIC, nsoil_water, &
    vgrid_soil_sb => new_vgrid('soil_layer_sb', ZAXIS_GENERIC, nsoil, &
                      longname='Soil (water) layers from SSE_ QUINCY', &
                      units='m', &
                      ! levels=  , &      ! is calculated from ubounds and dz
                      ! lbounds=  , &     ! is calculated from ubounds and dz
                      !ubounds=ubounds_soil_lay_water(:), & 
                      ubounds=soil_vgrid%ubounds(:), & 
                      !dz=dz_water(:))
                      dz=soil_vgrid%dz(:))
    CALL register_vgrid(vgrid_soil_sb)

  END SUBROUTINE Init_sse_config

#endif
END MODULE mo_sse_config_class
