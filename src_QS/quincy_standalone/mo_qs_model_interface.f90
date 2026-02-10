!> Interface for the QUINCY model to run one time step
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
MODULE mo_qs_model_interface

  USE mo_kind
  USE mo_exception

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: land_model_interface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qs_model_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> land model interface
  !!
  !! selection of the model via the namlist.slm - BASE_CTL:quincy_model_name
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE land_model_interface(tile                            , &
                                  buffer_nc_size                  , &
                                  buffer_nc_index                 , &
                                  l_forcing                       , &
                                  l_output                        , &
                                  l_output_as_ascii_text          , &
                                  l_output_interval_match         , &
                                  flux_trans_timestep_out_state   , &
                                  output_interval_pool_act        , &
                                  elapsed_seconds                 , &
                                  time_day, time_year)

    USE mo_jsb_class,                     ONLY: Get_model
    USE mo_jsb_tile_class,                ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                ONLY: t_jsb_process_task
    USE mo_qs_forcing,                    ONLY: forcing_options
    USE mo_jsb_model_class,               ONLY: t_jsb_model
    USE mo_quincy_model_config,           ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION    
    USE mo_jsb_process_class,             ONLY: A2L_, L2A_, SEB_, TURB_, SPQ_, HYDRO_, VEG_, HD_, Q_ASSIMI_, &
      &                                         Q_PHENO_, SB_, DFIRE_, PHYD_, Q_SYL_, Q_AGR_, &
      &                                         Get_process_name
    USE mo_hsm_class,                     ONLY: t_Message
    USE mo_qs_output_txt,                 ONLY: write_output_pool, write_output_flux 
    USE mo_qs_output_txt,                 ONLY: write_soil_forcing
    !> NOTE: writing of soil forcing still required if 'l_output_as_ascii_text = .FALSE.', i.e. if 'normal' output is written to netcdf
    USE mo_qs_output_netcdf4,             ONLY: write_netcdf4_quincy_buffered, accumulate_on_nc_acc_vars, &
      &                                         add_to_nc_buffer, reset_nc_acc_vars
    USE mo_qs_test_canopy_model,          ONLY: run_test_canopy_model
    USE mo_qs_test_radiation_model,       ONLY: run_test_radiation_model

    ! Process interfaces:
    USE mo_q_veg_turnover,                ONLY: update_veg_turnover
    USE mo_q_veg_dynamics,                ONLY: update_veg_dynamics
    USE mo_q_veg_growth,                  ONLY: update_veg_growth
    USE mo_q_veg_update_pools,            ONLY: update_veg_pools
    USE mo_q_veg_plant_characteristics,   ONLY: update_plant_characteristics
    USE mo_q_veg_products_decay,          ONLY: update_products_decay
    USE mo_q_veg_update_pools_on_harvest, ONLY: update_pools_on_harvest
    USE mo_q_assimi_interface,            ONLY: update_canopy_fluxes, update_ftranspiration_per_sl, &
      &                                         update_reset_q_assimi_fluxes, update_time_average_q_assimilation
    USE mo_spq_util,                      ONLY: reset_spq_fluxes
    USE mo_q_spq_update_turbulence,       ONLY: update_soil_turbulence
    USE mo_q_spq_update_physics,          ONLY: update_spq_physics
    USE mo_q_sb_rate_modifier,            ONLY: update_sb_rate_modifier
    USE mo_q_sb_ssm_main,                 ONLY: update_sb_simple_model
    USE mo_q_sb_jsm_main,                 ONLY: update_sb_jsm
    USE mo_q_sb_update_pools,             ONLY: update_sb_pools
    USE mo_sb_util,                       ONLY: reset_sb_fluxes, calculate_time_average_soilbiogeochemistry
    USE mo_q_pheno_process,               ONLY: update_q_phenology
    USE mo_phyd_process,                  ONLY: update_plant_hydraulics
    USE mo_q_syl_interface,               ONLY: update_sylviculture
    USE mo_q_agr_interface,               ONLY: update_cropland_process
    USE mo_dfire_interface,               ONLY: update_dfire_sf_fuel, update_dfire_sf_ignition, update_dfire_sf_spread, &
      &                                         update_dfire_sf_effects, update_dfire_sf_feedback
    USE mo_veg_util,                      ONLY: reset_veg_fluxes, calculate_time_average_vegetation, test_carbon_conservation
    USE mo_rad_interface,                 ONLY: update_surface_radiation, update_radiation_par, update_albedo, &
      &                                         update_time_average_radiation
    USE mo_turb_interface,                ONLY: update_humidity_scaling, update_roughness, aggregate_roughness
    USE mo_seb_interface,                 ONLY: update_surface_energy, update_asselin, &
      &                                         update_surface_fluxes, update_energy_balance
    USE mo_sse_interface,                 ONLY: update_soil_and_snow_temperature, update_soil_and_snow_properties
    USE mo_hydro_interface,               ONLY: update_surface_hydrology, update_soil_properties, &
      &                                         update_soil_hydrology, update_evaporation, &
      &                                         update_canopy_cond_unstressed, update_water_stress, update_canopy_cond_stressed, &
      &                                         update_snow_and_wet_fraction, update_water_balance

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(DFIRE_), DFIRE_MODEL_SF
    dsl4jsb_Use_config(PHYD_)
    dsl4jsb_Use_config(Q_SYL_)
    dsl4jsb_Use_config(Q_AGR_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile                           !< one tile with data structure for one lct
    INTEGER,                    INTENT(inout)     :: buffer_nc_index                !< current index in the nc output buffer
    INTEGER,                    INTENT(in)        :: buffer_nc_size                 !< size of the buffer used for nc output
    LOGICAL,                    INTENT(in)        :: l_forcing, &                   !< IF .TRUE. THEN update_forcing()
                                                     l_output, &                    !< TRUE if output is enabled at all (i.e., fluxes are accumulated)
                                                     l_output_as_ascii_text, &      !< TRUE if output to txt, else nc4
                                                     l_output_interval_match        !< TRUE if output may be created at the current timestep
    INTEGER,                    INTENT(in)        :: flux_trans_timestep_out_state  !< possible values: 0, 1, 2, 3; see quincy_qs_model_driver
    REAL(wp),                   INTENT(in)        :: output_interval_pool_act       !< output period in days; identical for pools and fluxes
    REAL(wp),                   INTENT(in)        :: elapsed_seconds                !< seconds elapsed since start of this simulation
    REAL(wp),                   INTENT(in)        :: time_day                       !< displays the 'real simulation time' in days; time_day == 1.0 after one day has passed
    INTEGER,                    INTENT(in)        :: time_year                      !< displays the 'real simulation time' in years; the 1st year of a simulation == 1 (not zero)

    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),          POINTER   :: model
    INTEGER                               :: ics, ice, iblk
    CHARACTER(len=*),           PARAMETER :: routine = TRIM(modname)//':land_model_interface'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(DFIRE_)
    dsl4jsb_Def_config(PHYD_)
    dsl4jsb_Def_config(Q_SYL_)
    dsl4jsb_Def_config(Q_AGR_)
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_config(PHYD_)
    dsl4jsb_Get_config(Q_SYL_)
    dsl4jsb_Get_config(Q_AGR_)

    ics  = model%options%ics
    ice  = model%options%ice
    iblk = model%options%iblk

    ! ------------------------------------------------------------------------------------------------------------
    ! Note: this is commented out because it would only be required in two cases that are currently not true
    !       1. it would be required if we would have a restart mechanism or
    !       2. if we initialise the variables belonging to the bgcms with values other than zero
    ! Furthermore: because some bgcm matrices are filled from forcing files in the SOIL only model,
    !              writing the variables to the matrices here would overwrite that information
    !!! store all bgcms into their matrices for this timestep
    !IF (ASSOCIATED(tile%bgcm_store)) THEN
    !  CALL tile%bgcm_store%Store_bgcms_in_matrices(tile%name, ics, ice, iblk)
    !ENDIF

    !> 2.0 Run the tasks defined for the specific SLM model
    !! 
    !! land plant soil canopy \n 
    !! test_canopy is a special case and runs for only one timestep
    !!
    SELECT CASE(model%config%qmodel_id)
    !! canopy, plant, or land ---------------------------------------------------------------------
    CASE(QCANOPY, QPLANT, QLAND)
      ! CANOPY, PLANT and LAND share many tasks
      ! QCANOPY: CANOPY model mode - basic set of tasks excluding the ones specific to PLANT or LAND
      !                              (calling only tasks to calculate biogeophysical processes)
      ! QPLANT:  PLANT  model mode - with additional tasks calculating vegetation biogeochemistry and dynamics
      ! QLAND:   LAND   model mode - adding tasks to calculate soil biogeochemistry
      !                              (on top of the tasks specific for the PLANT model mode)

      !! Assertions on combinations of processes on a QS tile
      ! In IQ Q_AGR_ will only run on crops, Q_SYL_ will only run on forests (usecase)
      ! ... in QS these should thus not be activated at the same tile
      IF (dsl4jsb_Config(Q_AGR_)%active .AND. dsl4jsb_Config(Q_SYL_)%active) THEN
        CALL finish( &
          & 'land_model_interface','Agriculture and sylviculture should not run ' &
          & // 'together for a tile. Please check.')
      END IF

      !! init/reset fluxes (with/to zero)
      !!
      CALL update_reset_q_assimi_fluxes(tile, model%options)
      IF (.NOT. model%config%use_soil_phys_jsbach) THEN
        CALL reset_spq_fluxes (tile, model%options)    ! mo_spq_util
      END IF
      CALL reset_veg_fluxes (tile, model%options)    ! mo_veg_util
      ! Soil biogeochemisty is only executed for LAND
      IF (model%config%qmodel_id == QLAND) THEN
        CALL reset_sb_fluxes (tile, model%options)   ! mo_sb_util; LAND only
      END IF

      !! tasks that might change area

      !... assert that syl is only active in land or plant mode
      IF (model%config%qmodel_id /= QLAND .AND. model%config%qmodel_id /= QPLANT) THEN
        IF (dsl4jsb_Config(Q_SYL_)%active) THEN
          PRINT *, 'Warning: sylviculture can currently only be active in land or plant mode, set active to .FALSE.'
          dsl4jsb_Config(Q_SYL_)%active = .FALSE.
        END IF
      END IF

      IF (tile%lcts(1)%lib_id > 0) THEN
        IF (dsl4jsb_Config(Q_AGR_)%active .AND. model%lctlib(tile%lcts(1)%lib_id)%CropFlag) THEN
           CALL update_cropland_process(tile, model%options)
        END IF
      END IF

      ! ... only executed for PLANT or LAND
      IF (model%config%qmodel_id == QPLANT .OR. model%config%qmodel_id == QLAND) THEN
        !!  LCC, LUC, disturbances
        !!
        IF (tile%lcts(1)%lib_id > 0) THEN
          IF (dsl4jsb_Config(Q_SYL_)%active .AND. dsl4jsb_Config(Q_SYL_)%harvest_active_in_qs) THEN
            IF (model%lctlib(tile%lcts(1)%lib_id)%ForestFlag) THEN
              CALL update_sylviculture(tile, model%options)        ! LAND or PLANT; if Q_SYL active
            END IF
          END IF
          IF ((dsl4jsb_Config(Q_AGR_)%active .AND. model%lctlib(tile%lcts(1)%lib_id)%CropFlag) &
              & .OR. (dsl4jsb_Config(Q_SYL_)%active .AND. model%lctlib(tile%lcts(1)%lib_id)%ForestFlag &
              &       .AND. dsl4jsb_Config(Q_SYL_)%harvest_active_in_qs)) THEN
            CALL update_pools_on_harvest(tile, model%options)      ! LAND or PLANT; if Q_SYL active
            CALL update_plant_characteristics(tile, model%options) ! LAND or PLANT; here only called if Q_SYL or Q_AGR are active
          END IF
        END IF

        !! disturbance
        !!   fire
        !!
        IF (dsl4jsb_Config(DFIRE_)%flag_dfire) THEN
          SELECT CASE (dsl4jsb_Config(DFIRE_)%dfire_model_id)
          CASE (DFIRE_MODEL_SF)
            CALL update_dfire_sf_fuel      (tile, model%options) ! PLANT or LAND (if DFIRE_)
            CALL update_dfire_sf_ignition  (tile, model%options) ! PLANT or LAND (if DFIRE_)
            CALL update_dfire_sf_spread    (tile, model%options) ! PLANT or LAND (if DFIRE_)
            CALL update_dfire_sf_effects   (tile, model%options) ! PLANT or LAND (if DFIRE_)
            CALL update_dfire_sf_feedback  (tile, model%options) ! PLANT or LAND (if DFIRE_)
          END SELECT
        END IF
      END IF

      !! tasks that do not change area
      !!
      ! Soil biogeochemisty is only executed for LAND
      IF(model%config%qmodel_id == QLAND) THEN
        CALL update_sb_rate_modifier(tile, model%options)        ! LAND ONLY
      ENDIF
      IF (model%config%use_soil_phys_jsbach) THEN
        CALL update_surface_radiation (tile, model%options)        ! jsbach RAD_
        CALL update_q_phenology (tile, model%options)
        CALL update_snow_and_wet_fraction(tile, model%options)     ! jsbach HYDRO_
        CALL update_soil_properties(tile, model%options)           ! jsbach HYDRO_
        CALL update_surface_energy(tile, model%options)            ! jsbach SEB_
        CALL update_evaporation(tile, model%options)               ! jsbach HYDRO_
        CALL update_surface_fluxes(tile, model%options)            ! jsbach SEB_
        CALL update_surface_hydrology(tile, model%options)         ! jsbach HYDRO_
        CALL update_snow_and_wet_fraction(tile, model%options)     ! jsbach HYDRO_
        CALL update_soil_and_snow_properties(tile, model%options)  ! jsbach SSE_
        CALL update_soil_and_snow_temperature(tile, model%options) ! jsbach SSE_
        CALL update_ftranspiration_per_sl(tile, model%options)     ! quincy Q_ASSIMI_
        CALL update_soil_hydrology(tile, model%options)            ! jsbach HYDRO_
        ! CALL update_canopy_cond_unstressed(tile, model%options)    ! jsbach HYDRO_ -- only used with jsbach in icon-land
        CALL update_water_stress(tile, model%options)              ! jsbach HYDRO_
        ! CALL update_canopy_cond_stressed(tile, model%options)      ! jsbach HYDRO_ -- only used with jsbach in icon-land
        ! -----------------------------------------
        ! only use with quincy in icon-land -------
        CALL update_radiation_par (tile, model%options)            ! jsbach RAD_
        CALL update_canopy_fluxes (tile, model%options)            ! quincy Q_ASSIMI_
        CALL update_humidity_scaling (tile, model%options)         ! jsbach TURB_
        ! -----------------------------------------
      END IF

      IF (.NOT. model%config%use_soil_phys_jsbach) THEN
        CALL update_humidity_scaling (tile, model%options)         ! jsbach TURB_
        CALL update_soil_turbulence (tile, model%options)          ! quincy SPQ_
        CALL update_surface_radiation (tile, model%options)        ! jsbach RAD_
        CALL update_q_phenology (tile, model%options)              ! quincy PHENO_
        CALL update_radiation_par (tile, model%options)            ! jsbach RAD_
        CALL update_canopy_fluxes (tile, model%options)            ! quincy Q_ASSIMI_        
        CALL update_ftranspiration_per_sl(tile, model%options)     ! quincy Q_ASSIMI_
        CALL update_spq_physics (tile, model%options)              ! quincy SPQ_
      END IF

      IF (model%config%qmodel_id == QPLANT .OR. model%config%qmodel_id == QLAND) THEN
        IF (dsl4jsb_Config(PHYD_)%use_plant_hydraulics) THEN
          CALL update_plant_hydraulics (tile, model%options)
        END IF
      END IF

      ! Vegetation dynamics and biogeochemisty are only calculated for PLANT or LAND
      IF (model%config%qmodel_id == QPLANT .OR. model%config%qmodel_id == QLAND) THEN
        IF(dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
          CALL update_products_decay(tile, model%options)        ! PLANT or LAND
        ENDIF
        CALL update_veg_turnover (tile, model%options)           ! PLANT or LAND
        CALL update_veg_dynamics (tile, model%options)           ! PLANT or LAND
        ! Soil biogeochemisty is only executed for LAND
        IF(model%config%qmodel_id == QLAND) THEN
          SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
          CASE("simple_1d")
            CALL update_sb_simple_model (tile, model%options)    ! LAND ONLY
          CASE('jsm')
            CALL update_sb_jsm (tile, model%options)             ! LAND ONLY
          END SELECT
        END IF
        CALL update_veg_growth (tile, model%options)             ! PLANT or LAND
      ELSE
        IF(dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
          CALL finish( &
            & 'land_model_interface','Using product pools only makes sense when running in PLANT or LAND mode. ' &
            & // ' Current mode: '//TRIM(model%config%quincy_model_name//'. Please check.'))
        END IF
      END IF

      CALL update_veg_pools (tile, model%options)
      CALL update_plant_characteristics(tile, model%options)

      ! Soil biogeochemisty is only executed for LAND
      IF(model%config%qmodel_id == QLAND) THEN
        CALL update_sb_pools (tile, model%options)               ! LAND ONLY
      END IF

      CALL update_roughness (tile, model%options)                ! jsbach (surface_roughness) TURB_
      CALL aggregate_roughness (tile, model%options)             ! jsbach (surface_roughness) TURB_
      CALL update_albedo(tile, model%options)                    ! jsbach + quincy RAD_

      IF (model%config%use_soil_phys_jsbach) THEN
        CALL update_water_balance(tile, model%options)             ! jsbach HYDRO_
        CALL update_asselin(tile, model%options)                   ! jsbach SEB_ | should always be the last task
        CALL update_energy_balance(tile, model%options)            ! jsbach SEB_ | is using 't' that is modified in 'update_asselin'
      END IF

      !! variable accumulation over time - daytime & moving average
      !!
      CALL update_time_average_radiation(tile, model%options)
      CALL update_time_average_q_assimilation(tile, model%options)
      CALL calculate_time_average_vegetation(tile, model%options)
      IF(model%config%qmodel_id == QLAND .OR. model%config%qmodel_id == QPLANT) THEN
        CALL calculate_time_average_soilbiogeochemistry(tile, model%options)
      END IF

      !! Carbon conservation test
      IF(model%config%qmodel_id == QLAND) THEN
        CALL test_carbon_conservation(tile, model%options%dtime)
      END IF

      ! Soil forcing can only be written for LAND
      IF(model%config%qmodel_id == QLAND) THEN
        !! create soil forcing for quincy-jsm standalone
        !! the output-routine is called each timestep, hence, output interval is timestep
        !! output litterfall and soilphysics from a static simulation
        IF (model%config%flag_generate_soil_forcing_stat) THEN
          IF (.NOT. model%config%l_transient_spinup .AND. &
              & time_year >= model%config%output_start_first_day_year &
              & .AND. time_year <= model%config%output_end_last_day_year) THEN
            CALL write_soil_forcing(tile, model%options, time_day)
          END IF
        END IF
        !! output litterfall and soilphysics from a transient fluxnet/3 diagnostics simulation
        IF (model%config%flag_generate_soil_forcing_trans .AND. l_output) THEN
          ! write soil forcing only for the last full forcing period of the particular fluxnet site (where static forcing is available)
          IF (flux_trans_timestep_out_state == 2) THEN
            CALL write_soil_forcing(tile, model%options, time_day)
          END IF
        END IF
        !! output litterfall and soilphysics from a transient cruncep-type site-set simulation (for years 1901 to ...)
        IF (model%config%flag_generate_soil_forcing_trans .AND. l_output) THEN
          ! write only after spinup and if not a transient fluxnet/3 diagnostics simulation
          IF (.NOT. model%config%l_transient_spinup &
              & .AND. .NOT. model%config%fluxnet_type_transient_timestep_output) THEN
            CALL write_soil_forcing(tile, model%options, time_day)
          END IF
        END IF
      END IF

    !! soil ---------------------------------------------------------------------------------------
    CASE(QSOIL)

      !! init/reset fluxes (with/to zero)
      !!
      CALL reset_sb_fluxes            (tile, model%options)   ! mo_sb_util
      CALL reset_spq_fluxes           (tile, model%options)   ! mo_spq_util

      !! tasks that do not change area
      !!
      CALL update_sb_rate_modifier    (tile, model%options)
      SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
      CASE("simple_1d")
        CALL update_sb_simple_model (tile, model%options)
      CASE('jsm')
        CALL update_sb_jsm          (tile, model%options)
      ENDSELECT
      CALL update_sb_pools            (tile, model%options)

      !! variable accumulation over time - daytime & moving average
      !!
      CALL update_time_average_radiation(tile, model%options)
      CALL update_time_average_q_assimilation(tile, model%options)
      CALL calculate_time_average_vegetation(tile, model%options)
      CALL calculate_time_average_soilbiogeochemistry(tile, model%options)

    CASE(Q_TEST_CANOPY)
    !! test_canopy --------------------------------------------------------------------------------
      CALL run_test_canopy_model  (tile, model%options)

    CASE(Q_TEST_RADIATION)
    !! test_radiation -----------------------------------------------------------------------------
      CALL run_test_radiation_model (tile, model%options)

    CASE DEFAULT
      CALL finish( &
        & 'land_model_interface','Invalid SLM model name. Try either land, plant, soil, canopy, test_canopy or test_radiation.')
    END SELECT  ! select SLM model !

    !! push matrices into bgcm variables
    IF (ASSOCIATED(tile%bgcm_store)) THEN
      CALL tile%bgcm_store%Write_stored_matrices_to_bgcms(tile%name, ics, ice, iblk)
    ENDIF

    !>  2.1 Run output routines for land, plant, soil, canopy
    !>
    !>    the test_canopy and test_radiation model modes feature their own output statements
    !>
    IF (l_output) THEN 
      SELECT CASE(model%config%qmodel_id)
      CASE(QLAND, QPLANT, QSOIL, QCANOPY)

        !! output pool: write output in certain intervals, see mo_qs_output_txt
        !!  note: weekly interval ends with 364 in a year, and starts with 7 in the next year -> 8 days
        !!  note: monthly interval is defined as 30 days, thus, ends with 360 and start in the next year with 30 -> 35 days
        !! 
        IF (l_output_as_ascii_text) THEN
          IF (l_output_interval_match) THEN
            CALL write_output_pool(tile, model%options, time_day)
          ENDIF
          !! output flux: write output in certain intervals ('output_interval_pool_act'),
          !!              see mo_qs_output_txt
          !!  if output is enabled, flux variables are accumulated
          !!  output is written when the interval is a match with current timestep
          !! 
          CALL write_output_flux(tile, model%options, time_day, output_interval_pool_act, &
            &                    write_flux_txt_output=l_output_interval_match, &
            &                    set_zero_flux_acc_variables=.FALSE.)
        ELSE
          CALL accumulate_on_nc_acc_vars(tile, model%options, output_interval_pool_act)
          ! write to nc4 files
          IF (l_output_interval_match) THEN
            CALL add_to_nc_buffer(tile, model%options, buffer_nc_index, elapsed_seconds)
            CALL reset_nc_acc_vars(tile, model%options)
            buffer_nc_index = buffer_nc_index + 1
          ENDIF
          
          ! If we reach the end of the buffer -> write the buffer to the netcdf file
          IF (buffer_nc_index  == buffer_nc_size  + 1) THEN
            CALL write_netcdf4_quincy_buffered(tile, model%options, buffer_nc_size)
            ! Reset the buffer index after writing
            buffer_nc_index = 1
          ENDIF 
        ENDIF
      END SELECT  ! select SLM model !
    ENDIF

  END SUBROUTINE land_model_interface

END MODULE mo_qs_model_interface
