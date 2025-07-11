!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!> 
!> @brief Initialization of the the dfire memory.
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
!> @par Revision History
!>
MODULE mo_dfire_init
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message
  USE mo_jsb_control,             ONLY: debug_on

  USE mo_dfire_constants,         ONLY: df_sf_fire_duration_thd
  USE mo_dfire_constants,         ONLY: ratio_as_dead_fuel
#ifndef __QUINCY_STANDALONE__
  USE mo_jsb_time,                ONLY: timestep_in_days
#endif
  USE mo_jsb_process_class,     ONLY: DFIRE_

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dfire_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_init'

  TYPE :: t_dfire_init_vars
    REAL(wp), POINTER ::      &
        pop_dens( :,:) => NULL(), &
        lightning_events(:,:,:) => NULL()
  END TYPE t_dfire_init_vars

  TYPE(t_dfire_init_vars) :: dfire_init_vars

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Intialize disturbance fire process (after memory has been set up)
  !>
  SUBROUTINE dfire_init(tile)
    USE mo_jsb_class,               ONLY: get_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,         ONLY: t_jsb_model

    USE mo_dfire_constants,         ONLY: DFIRE_PRECIP_METHOD_DAYTIME, & !< Method to compute precipitation
                                          DFIRE_PRECIP_METHOD_APROX_MAVG, &
                                          DFIRE_PRECIP_METHOD_MAVG

    dsl4jsb_Use_config(DFIRE_)
    dsl4jsb_Use_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER              :: model
    !INTEGER                                      :: iblk, ics, ice, nc, npft
    REAL(wp)                                     :: dtime
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':dfire_init'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(DFIRE_)
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: dfire_pm_id                                  !< Selected precipitation method
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real2D_onDomain      :: mask_is_fire_intense_enough 
    dsl4jsb_Real2D_onDomain      :: precip_daytime
    dsl4jsb_Real2D_onDomain      :: precip_qc_mavg
    dsl4jsb_Real2D_onDomain      :: precip_daily_mmd
    dsl4jsb_Real2D_onDomain      :: precip_mavg

    dsl4jsb_Real2D_onDomain      :: t_air_min 
    dsl4jsb_Real2D_onDomain      :: t_air_max

    dsl4jsb_Real2D_onDomain      :: nesterov_index
    dsl4jsb_Real2D_onDomain      :: fdi
    dsl4jsb_Real2D_onDomain      :: ros_f
    dsl4jsb_Real2D_onDomain      :: ros_b
    dsl4jsb_Real2D_onDomain      :: mean_fire_area
    dsl4jsb_Real2D_onDomain      :: fire_duration 

    dsl4jsb_Real2D_onDomain      :: area_burned
    dsl4jsb_Real2D_onDomain      :: area_burned_fract
    dsl4jsb_Real3D_onDomain      :: dead_frac_consumption
    dsl4jsb_Real2D_onDomain      :: i_surface

    dsl4jsb_Real3D_onDomain      :: fuel_moisture

    dsl4jsb_Real2D_onDomain      :: e_n_ignitions
    dsl4jsb_Real2D_onDomain      :: e_n_humans
    dsl4jsb_Real2D_onDomain      :: e_n_lightning
    dsl4jsb_Real2D_onDomain      :: ab_moisture_diag

    dsl4jsb_Real2D_onDomain      :: reaction_intensity
    dsl4jsb_Real2D_onDomain      :: prop_flux_ratio
    dsl4jsb_Real2D_onDomain      :: wind_factor
    dsl4jsb_Real2D_onDomain      :: eff_heating_number
    dsl4jsb_Real2D_onDomain      :: preig_heat

    dsl4jsb_Real2D_onDomain      :: opt_react_velo
    dsl4jsb_Real2D_onDomain      :: net_fuel_load
    dsl4jsb_Real2D_onDomain      :: moist_damp_coeff
    dsl4jsb_Real2D_onDomain      :: mineral_damp_coeff

    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_1hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_10hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_100hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_1000hr
    dsl4jsb_Real3D_onDomain      :: live_fuel_consumption_1hr

    dsl4jsb_Real2D_onDomain      :: distance_travelled
    dsl4jsb_Real2D_onDomain      :: l2b_ratio

    dsl4jsb_Real3D_onDomain      :: frac_consump_fcl
    dsl4jsb_Real3D_onDomain      :: ic_pop_dens_monthly
    dsl4jsb_Real3D_onDomain      :: ic_lightning_events_monthly
    ! ----------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    !iblk    = options%iblk
    !ics     = options%ics
    !ice     = options%ice
    !nc      = options%nc
    !dtime   = options%dtime
    dtime=timestep_in_days(model_id=tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(DFIRE_)) RETURN
    !IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => get_model(tile%owner_model_id)
    !hgrid  => Get_grid(model%grid_id)
    !nblks  =  hgrid%nblks
    !nproma =  hgrid%nproma
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onDomain(DFIRE_, nesterov_index)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, fdi)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, ros_f)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, ros_b)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, mean_fire_area)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, fire_duration)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, e_n_ignitions)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, e_n_humans)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, e_n_lightning)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, ab_moisture_diag)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, reaction_intensity)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, prop_flux_ratio)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, wind_factor)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, eff_heating_number)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, preig_heat)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, opt_react_velo)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, net_fuel_load)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, moist_damp_coeff)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, mineral_damp_coeff)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, mask_is_fire_intense_enough)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, precip_qc_mavg)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, precip_daily_mmd)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, precip_mavg)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, precip_daytime)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, t_air_min)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, t_air_max)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, area_burned)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, area_burned_fract)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_frac_consumption)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, i_surface)

    dsl4jsb_Get_var3D_onDomain(DFIRE_, fuel_moisture)

    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_fuel_consumption_1hr)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_fuel_consumption_10hr)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_fuel_consumption_100hr)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_fuel_consumption_1000hr)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, live_fuel_consumption_1hr)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, distance_travelled)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, l2b_ratio)

    dsl4jsb_Get_var3D_onDomain(DFIRE_, frac_consump_fcl)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, ic_pop_dens_monthly)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, ic_lightning_events_monthly)
    ! ----------------------------------------------------------------------------------------------------- !
    nesterov_index = 0._wp
    fdi = 0._wp
    ros_f = 0._wp
    ros_b = 0._wp
    mean_fire_area = 0._wp
    fire_duration = 0._wp

    e_n_ignitions = 0._wp
    e_n_humans = 0._wp
    e_n_lightning = 0._wp
    ab_moisture_diag = 0._wp

    mask_is_fire_intense_enough = 0._wp !< TRUE=1._wp, FALSE=0._wp
    precip_daytime = 0._wp
    precip_qc_mavg = 0._wp
    precip_daily_mmd = 0._wp
    precip_mavg = 0._wp

    t_air_min = 0._wp
    t_air_max = 0._wp

    area_burned = 0._wp
    area_burned_fract = 0._wp
    dead_frac_consumption = 0._wp
    i_surface = 0._wp
    dead_fuel_consumption_1hr = 0._wp
    dead_fuel_consumption_10hr = 0._wp
    dead_fuel_consumption_100hr = 0._wp
    dead_fuel_consumption_1000hr = 0._wp
    live_fuel_consumption_1hr = 0._wp

    fuel_moisture = 0._wp

    reaction_intensity = 0._wp
    prop_flux_ratio = 0._wp
    wind_factor = 0._wp
    eff_heating_number = 0._wp
    preig_heat = 0._wp

    distance_travelled = 0._wp
    l2b_ratio = 0._wp
    frac_consump_fcl = 0._wp

    df_sf_fire_duration_thd = dtime / 60._wp ! convert seconds to minutes

    !< Chose daily precipitation. Those are calculated in 3 different ways.
    dfire_pm_id = dsl4jsb_Config(DFIRE_)%dfire_precip_method_id
    IF (dfire_pm_id == DFIRE_PRECIP_METHOD_DAYTIME) THEN
      !< daytime
      dsl4jsb_var_ptr(DFIRE_,precip_daily) => dsl4jsb_var_ptr(DFIRE_,precip_daytime)
    ELSE IF (dfire_pm_id == DFIRE_PRECIP_METHOD_MAVG) THEN
      !< moving average
      dsl4jsb_var_ptr(DFIRE_,precip_daily) => dsl4jsb_var_ptr(DFIRE_,precip_mavg)
    ELSE IF (dfire_pm_id == DFIRE_PRECIP_METHOD_APROX_MAVG) THEN
      !< aproximation to moving average (quincy default method)
      dsl4jsb_var_ptr(DFIRE_,precip_daily) => dsl4jsb_var_ptr(DFIRE_,precip_qc_mavg)
    ENDIF

    !< Initialize ratio of sb pools to dead fuel
    ratio_as_dead_fuel % polymeric      = dsl4jsb_Config(DFIRE_)%ratio_polymeric_litter_to_dead_fuel
    ratio_as_dead_fuel % dom_assoc      = dsl4jsb_Config(DFIRE_)%ratio_dom_assoc_to_dead_fuel
    ratio_as_dead_fuel % dom            = dsl4jsb_Config(DFIRE_)%ratio_dom_to_dead_fuel
    ratio_as_dead_fuel % soluable       = dsl4jsb_Config(DFIRE_)%ratio_soluable_litter_to_dead_fuel
    ratio_as_dead_fuel % residue        = dsl4jsb_Config(DFIRE_)%ratio_residue_to_dead_fuel
    ratio_as_dead_fuel % residue_assoc  = dsl4jsb_Config(DFIRE_)%ratio_residue_assoc_to_dead_fuel

#ifdef __QUINCY_STANDALONE__
    ! quincy only
    CALL dfire_read_init_vars(tile)
#else
    ! IQ only
    ! call only at root tile, i.e., the only tile without associated parent tile (to avoid unnessary i/o)
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL dfire_read_init_vars(tile)
    END IF

    !CALL spq_init_ic_bc(tile)
    ic_pop_dens_monthly = dfire_init_vars%pop_dens
    ic_lightning_events_monthly = dfire_init_vars%lightning_events

    IF (tile%Is_last_process_tile(DFIRE_)) THEN
    !  CALL spq_finalize_init_vars()
    END IF
#endif

  END SUBROUTINE dfire_init


#ifdef __QUINCY_STANDALONE__
  SUBROUTINE dfire_read_init_vars(tile)
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(DFIRE_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(DFIRE_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = modname//':dfire_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting DFIRE_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels

    ALLOCATE( &
      & dfire_init_vars%pop_dens               (nproma, nblks),       &
      & dfire_init_vars%lightning_events       (nproma, nblks, 12),   &
      & )

    dfire_init_vars%pop_dens              (:,:) = missval
    dfire_init_vars%lightning_events    (:,:,:) = missval

    dsl4jsb_Get_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !

    dfire_init_vars%pop_dens(:,:) = dsl4jsb_Config(DFIRE_)%pop_dens
    dfire_init_vars%lightning_events(:,:,:) = dsl4jsb_Config(DFIRE_)%lightning_events
  END SUBROUTINE dfire_read_init_vars
#else
  ! ======================================================================================================= !
  !> read DFIRE_ init vars from input file
  !>
  SUBROUTINE dfire_read_init_vars(tile)
    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(DFIRE_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(DFIRE_)

    REAL(wp), POINTER           :: ptr_2D(:, :)
    REAL(wp), POINTER           :: ptr_2D_time(:, :, :)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb
    TYPE(t_input_file)          :: input_file
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = modname//':dfire_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting DFIRE_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()

    ALLOCATE( &
        & dfire_init_vars%pop_dens               (nproma, nblks),       &
        & dfire_init_vars%lightning_events       (nproma, nblks, 12)    & ! months in 1 year
        & )

    dfire_init_vars%pop_dens              (:,:) = missval
    dfire_init_vars%lightning_events    (:,:,:) = missval

    dsl4jsb_Get_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !

    ! Read data from populations density file
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DFIRE_)%ic_popdens_filename), model%grid_id)

    !> pop dens
    ptr_2D => input_file%Read_2d(  &
        & variable_name='Band1', &
        & fill_array = dfire_init_vars%pop_dens)
    ! Clean: set negative values to 0
    dfire_init_vars%pop_dens = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    CALL input_file%Close()


    ! Read data from sso file
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DFIRE_)%ic_lightning_events_filename), model%grid_id)

    !>
    ptr_2D_time => input_file%Read_2D_time(  &
        & variable_name='HRMC_COM_FR', &
        & fill_array = dfire_init_vars%lightning_events, &
        & start_time_step = 1, end_time_step = 12 )
    ! Clean: set negative values to 0
    dfire_init_vars%lightning_events = MERGE(ptr_2D_time, 0._wp, ptr_2D_time >= 0._wp)

    CALL input_file%Close()

  END SUBROUTINE dfire_read_init_vars
#endif

#endif
END MODULE mo_dfire_init
