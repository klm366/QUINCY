!> Summary: interface to the disturb process (disturbances: currently fire and wind)
!>
!>### Contains the interface to the disturb process (disturbances: currently fire and wind)
!>
!> Copyright (c) 2017 Rainer Schneck
!> Copyright (c) 2018-2021 Julia Nabel
!>
!>#### Revision History
!>
!> First version                                by Rainer Schneck (2017-05-04)
!> C conservation test                          by Julia Nabel (2018)
!> Disassembling of tasks and ongoing changes   by Julia Nabel (2021-10-29)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
!>
MODULE mo_disturb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_carbon_constants,   ONLY: molarMassCO2_kg, sec_per_day

  ! Use of processes in this module (Get_disturb_memory and Get_disturb_config)
  dsl4jsb_Use_processes DISTURB_, FUEL_, ASSIMI_, A2L_, CARBON_, SEB_, FLCC_, WLCC_
  dsl4jsb_Use_processes PHENO_, HYDRO_    ! <Rackhun Son>

  ! Use of process configurations (t_disturb_config)
  dsl4jsb_Use_config(DISTURB_)
  dsl4jsb_Use_config(FLCC_)
  dsl4jsb_Use_config(WLCC_)
  dsl4jsb_Use_config(CARBON_)

  ! Use of process memories (t_disturb_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(DISTURB_)
  dsl4jsb_Use_memory(FUEL_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(PHENO_)    ! <Rackhun Son>
  dsl4jsb_Use_memory(HYDRO_)    ! <Rackhun Son>

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_disturb_tasks !,t_disturb_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_interface'

  !> Type definition for damaged area task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_damaged_area
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_damaged_area    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_damaged_area !< Aggregates computed task variables
  END TYPE tsk_damaged_area

  !> Constructor interface for damaged area task
  INTERFACE tsk_damaged_area
    PROCEDURE Create_task_damaged_area     !< Constructor function for task
  END INTERFACE tsk_damaged_area

  !> Type definition for natural_disturbances task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_natural_disturbances
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_natural_disturbances    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_natural_disturbances !< Aggregates computed task variables
  END TYPE tsk_natural_disturbances

  !> Constructor interface for natural_disturbances task
  INTERFACE tsk_natural_disturbances
    PROCEDURE Create_task_natural_disturbances     !< Constructor function for task
  END INTERFACE tsk_natural_disturbances

CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for damaged area task
  !
  FUNCTION Create_task_damaged_area(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_damaged_area::return_ptr)
    CALL return_ptr%Construct(name='damaged_area', process_id=DISTURB_, owner_model_id=model_id)

  END FUNCTION Create_task_damaged_area

  ! ====================================================================================================== !
  !
  !> Constructor for C_natural_disturbance task
  !
  FUNCTION Create_task_natural_disturbances(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr
    ! -------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_natural_disturbances::return_ptr)
    CALL return_ptr%Construct(name='natural_disturbances', process_id=DISTURB_, owner_model_id=model_id)

  END FUNCTION Create_task_natural_disturbances

  ! ====================================================================================================== !
  !
  !> Register tasks for disturb process
  !
  SUBROUTINE Register_disturb_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id
    ! -------------------------------------------------------------------------------------------------- !
    CALL this%Register_task(tsk_damaged_area(model_id))
    CALL this%Register_task(tsk_natural_disturbances(model_id))

  END SUBROUTINE Register_disturb_tasks


  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "damaged_area"
  !>
  !> Task "damaged_area" calculates the damage of vegetation due to fire and windbreak
  !
  SUBROUTINE update_damaged_area(tile, options)
    USE mo_disturb_process,       ONLY: burned_fract_jsbach, broken_woody_fract_jsbach, get_relative_humidity_air
    USE mo_disturb_process,       ONLY: burned_fract_dl ! <Rackhun Son>  
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_time,              ONLY: is_time_experiment_start, timesteps_per_day, is_newday
    USE mo_disturbance_constants, ONLY: persist_rel_hum, persist_wind_10m

    USE mo_jsb_time_iface,        ONLY: get_year_day    ! <Rackhun Son>
    USE mo_time_config,           ONLY: time_config     ! <Rackhun Son>
    !!!! These will be removed after the FUEL issue
    USE mo_fuel_process,          ONLY: calc_fuel_jsbach
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    INTEGER :: iblk, nc, ics, ice
    INTEGER :: iyr, imo, idy, day_of_year, diff_year, i          ! <Rackhun Son>
    REAL(wp) :: clim_ratio, yearly_ratio                         ! average ratio for climatology <Rackhun Son>

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_damaged_area'

    ! Declare process memories
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_memory(FUEL_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(PHENO_)                       ! <Rackhun Son>
    dsl4jsb_Def_memory(HYDRO_)                       ! <Rackhun Son>
    dsl4jsb_Def_config(DISTURB_)
    dsl4jsb_Def_config(CARBON_)

    dsl4jsb_Real2D_onChunk :: fuel                   ! Amount of fuel on which the fire was based.
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: wind_10m               ! 10m wind speed [m/s]
    dsl4jsb_Real2D_onChunk :: max_wind_10m_act
    dsl4jsb_Real2D_onChunk :: prev_day_max_wind_10m
    dsl4jsb_Real2D_onChunk :: max_wind_10m
    dsl4jsb_Real2D_onChunk :: q_rel_air_climbuf
    dsl4jsb_Real2D_onChunk :: q_rel_air_climbuf_yDay

    !!!! These will be removed after the FUEL issue
    dsl4jsb_Real2D_onChunk :: c_acid_ag1_ta
    dsl4jsb_Real2D_onChunk :: c_water_ag1_ta
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag1_ta
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag1_ta
    dsl4jsb_Real2D_onChunk :: c_acid_ag2_ta
    dsl4jsb_Real2D_onChunk :: c_water_ag2_ta
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag2_ta
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag2_ta
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! axuxiliary input variables for fire
    dsl4jsb_Real2D_onChunk :: t_air                  ! temperature [K]   <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: rain                   ! precipitation <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: lai                    ! LAI <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: fract_fpc_max          ! Maximum fraction of tile that can be covered by leave <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: pfts_fract             ! fraction of PFTs <Rackhun Son>
    !dsl4jsb_Real2D_onChunk :: fract_snow_soil        ! snow fraction on soil <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: fract_snow             ! Snow fraction of tile <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: elevation              ! elevation <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl            ! water content in soil layers [m] <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: soil_depth_sl          ! soil depth within each soil layer [m] <Rackhun Son>

    dsl4jsb_Real2D_onChunk :: t_air_ave              ! daily average of temperature <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: q_ave                  ! daily average of specific humidity <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: q_rel_ave              ! daily average of relative humidity <Rackhun Son>
    !dsl4jsb_Real2D_onChunk :: ws10_ave               ! daily average of wind speed <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl1_ave               ! daily average of volumetric soil water content in the 1st layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl2_ave               ! daily average of volumetric soil water content in the 2nd layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl3_ave               ! daily average of volumetric soil water content in the 3rd layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl4_ave               ! daily average of volumetric soil water content in the 4th layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: lai_ave                ! daily average of LAI <Rackhun Son>

    dsl4jsb_Real3D_onChunk :: t_air_clim             ! daily climatology of temperature <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: q_clim                 ! daily climatology of specific humidity <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: q_rel_clim             ! daily climatology of relative humidity <Rackhun Son>
    !dsl4jsb_Real3D_onChunk :: ws10_clim              ! daily climatology of wind speed <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl1_clim              ! daily climatology of volumetric soil water content in the 1st layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl2_clim              ! daily climatology of volumetric soil water content in the 2nd layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl3_clim              ! daily climatology of volumetric soil water content in the 3rd layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl4_clim              ! daily climatology of volumetric soil water content in the 4th layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: lai_clim               ! daily climatology of LAI <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: fuel_yearly            ! yearly mean of Fuel <Rackhun Son>

    dsl4jsb_Real3D_onChunk :: climate_lstm_hidden    ! LSTM hidden/cell <Rackhun>
    dsl4jsb_Real3D_onChunk :: climate_lstm_cell      ! LSTM hidden/cell <Rackhun>
    dsl4jsb_Real3D_onChunk :: vegetation_lstm_hidden ! LSTM hidden/cell <Rackhun>
    dsl4jsb_Real3D_onChunk :: vegetation_lstm_cell   ! LSTM hidden/cell <Rackhun>

    dsl4jsb_Real2D_onChunk :: burned_fract
    dsl4jsb_Real2D_onChunk :: damaged_fract

    dsl4jsb_Real2D_onChunk :: co2flux_fire_all_2_atm
    dsl4jsb_Real2D_onChunk :: cflux_dist_green_2_soil
    dsl4jsb_Real2D_onChunk :: cflux_dist_woods_2_soil

    dsl4jsb_Real2D_onChunk :: cover_fract_pot

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & q_rel_air,                     &   
      & swl1,                          &   ! volumetric soil water content in the 1st layer <Rackhun Son>
      & swl2,                          &   ! volumetric soil water content in the 2nd layer <Rackhun Son>
      & swl3,                          &   ! volumetric soil water content in the 3rd layer <Rackhun Son>
      & swl4                               ! volumetric soil water content in the 4th layer <Rackhun Son>                

    REAL(wp), DIMENSION(options%nc, 9) :: &
      & modified_pfts_fract                ! modified pft fractions for Fire <Rackhun Son>

    ! Declare pointers to variables in memory <Rackhun Son>
    REAL(wp), POINTER      :: lightning_climatology(:)
    REAL(wp), POINTER      :: population_density(:)
    REAL(wp), POINTER      :: t_air_local_clim(:)
    REAL(wp), POINTER      :: q_rel_local_clim(:)
    REAL(wp), POINTER      :: q_local_clim(:)
    REAL(wp), POINTER      :: ws10_local_clim(:)
    REAL(wp), POINTER      :: swvl1_local_clim(:)
    REAL(wp), POINTER      :: swvl2_local_clim(:)
    REAL(wp), POINTER      :: swvl3_local_clim(:)
    REAL(wp), POINTER      :: swvl4_local_clim(:)
    REAL(wp), POINTER      :: lai_local_clim(:)

    !LOGICAL  :: pft_mask ! True for pfts to which the disturbance should be applied.
    REAL(wp) :: persist ! Factor (as fraction) determining the persistance of q_rel_air_climbuf
                        ! relative to the time-step-actual q_rel_air
    REAL(wp) :: fire_minimum, fire_tau
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    iyr     = options%current_datetime%date%year   ! <Rackhun Son> 
    imo     = INT(options%current_datetime%date%month)   ! <Rackhun Son> 
    idy     = options%current_datetime%date%day     ! <Rackhun Son> 

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(DISTURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process configurations
    dsl4jsb_Get_config(DISTURB_)
    dsl4jsb_Get_config(CARBON_)

    ! Get process memories
    dsl4jsb_Get_memory(DISTURB_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(FUEL_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(PHENO_)    ! <Rackhun Son>
    dsl4jsb_Get_memory(HYDRO_)    ! <Rackhun Son>

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(FUEL_,      fuel)                     ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,       q_air)                    ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,       press_srf)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,       wind_10m)                 ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,       t_unfilt)                 ! in
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   max_wind_10m_act)         ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   prev_day_max_wind_10m)    ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   max_wind_10m)             ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_rel_air_climbuf)        ! inout
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_rel_air_climbuf_yDay)   ! out

    dsl4jsb_Get_var2D_onChunk(A2L_,       t_air)                    ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(A2L_,       rain)                     ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(PHENO_,     lai)                      ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(PHENO_,     fract_fpc_max)            ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(HYDRO_,     pfts_fract)               ! in <Rackhun Son>
    !dsl4jsb_Get_var2D_onChunk(HYDRO_,     fract_snow_soil)         ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(HYDRO_,     fract_snow)               ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(HYDRO_,     elevation)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(HYDRO_,     wtr_soil_sl)              ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(HYDRO_,     soil_depth_sl)            ! in <Rackhun Son>
    
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   t_air_ave)                ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_ave)                    ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_rel_ave)                ! in <Rackhun Son>
    !dsl4jsb_Get_var2D_onChunk(DISTURB_,   ws10_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl1_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl2_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl3_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl4_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   lai_ave)                  ! in <Rackhun Son>

    dsl4jsb_Get_var3D_onChunk(DISTURB_,   t_air_clim)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   q_clim)                   ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   q_rel_clim)               ! in <Rackhun Son>
    !dsl4jsb_Get_var3D_onChunk(DISTURB_,   ws10_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl1_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl2_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl3_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl4_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   lai_clim)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   fuel_yearly)              ! in <Rackhun Son>

    dsl4jsb_Get_var3D_onChunk(DISTURB_,   climate_lstm_hidden)      ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   climate_lstm_cell)        ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   vegetation_lstm_hidden)   ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   vegetation_lstm_cell)     ! in <Rackhun Son>
    
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   burned_fract)             ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   damaged_fract)            ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,    co2flux_fire_all_2_atm)   ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_green_2_soil)  ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_woods_2_soil)  ! out
    
    !!!! These will be removed after the FUEL issue 
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_ag1_ta)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_ag1_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_ag1_ta)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_ag1_ta)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_ag2_ta)            ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_ag2_ta)           ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_ag2_ta)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_ag2_ta)      ! in
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    dsl4jsb_Get_var2D_onChunk(ASSIMI_,    cover_fract_pot)          ! in
    ! -------------------------------------------------------------------------------------------------- !

    q_rel_air(:) = get_relative_humidity_air(nc, q_air(:), t_unfilt(:), press_srf(:))

    ! For each time step update the climate buffer

    ! R: This may be later in a climbuf process for dynveg and disturbance.
    ! Maximum daily wind speed and relative air humidity smoothed in time
    persist = persist_rel_hum ** (1._wp/REAL(timesteps_per_day(options%dtime),wp))
    q_rel_air_climbuf(:) = q_rel_air_climbuf(:) * persist  +  MIN(q_rel_air(:),100._wp) * (1._wp - persist)

    IF (.NOT. is_newday(options%current_datetime, options%dtime) .OR. is_time_experiment_start(options%current_datetime)) THEN
      WHERE (wind_10m(:) > max_wind_10m_act(:))
         max_wind_10m_act(:) = wind_10m(:)
      END WHERE
    ELSE
      prev_day_max_wind_10m(:) = max_wind_10m_act(:)
      max_wind_10m(:) = max_wind_10m(:) * persist_wind_10m + max_wind_10m_act(:) * (1._wp - persist_wind_10m)
      max_wind_10m_act(:) = wind_10m(:)
      q_rel_air_climbuf_yDay(:) = q_rel_air_climbuf(:)
    END IF

    ! calculate volumetric soil water content <Rackhun Son> 
    !swl1(:) = wtr_soil_sl(:,1) / (soil_depth_sl(:,1)+1e-8_wp)
    !swl2(:) = wtr_soil_sl(:,2) / (soil_depth_sl(:,2)+1e-8_wp)
    !swl3(:) = wtr_soil_sl(:,3) / (soil_depth_sl(:,3)+1e-8_wp)
    !swl4(:) = wtr_soil_sl(:,4) / (soil_depth_sl(:,4)+1e-8_wp)

    DO i = 1, size(soil_depth_sl, dim=1)
      IF (soil_depth_sl(i,1) > 0._wp) THEN
        swl1(i) = wtr_soil_sl(i,1) / soil_depth_sl(i,1)
      ELSE
        swl1(i) = 0._wp
      ENDIF

      IF (soil_depth_sl(i,2) > 0._wp) THEN
        swl2(i) = wtr_soil_sl(i,2) / soil_depth_sl(i,2)
      ELSE
        swl2(i) = 0._wp
      ENDIF

      IF (soil_depth_sl(i,3) > 0._wp) THEN
        swl3(i) = wtr_soil_sl(i,3) / soil_depth_sl(i,3)
      ELSE
        swl3(i) = 0._wp
      ENDIF

      IF (soil_depth_sl(i,4) > 0._wp) THEN
        swl4(i) = wtr_soil_sl(i,4) / soil_depth_sl(i,4)
      ELSE
        swl4(i) = 0._wp
      ENDIF
    END DO

    ! calculate daily statistics <Rackhun Son>
    ! get day of year for cliamtology calculation
    day_of_year = INT(get_year_day(options%current_datetime))

    ! skip 2/29 (2/29 is 60th at proleptic gregorian)
    IF (mod(iyr,4) /= 0) THEN
      IF (day_of_year >= 60) day_of_year = day_of_year+1
    ENDIF

    !where (q_rel_air > 100) q_rel_air = 100.

    t_air_ave(:) = t_air_ave(:) + t_air(:)    *(options%dtime/60._wp/60._wp/24._wp)
    q_ave(:)     = q_ave(:)     + q_air(:)    *(options%dtime/60._wp/60._wp/24._wp)
    q_rel_ave(:) = q_rel_ave(:) + q_rel_air(:)*(options%dtime/60._wp/60._wp/24._wp)
    !ws10_ave(:)  = ws10_ave(:)  + wind_10m(:) *(options%dtime/60._wp/60._wp/24._wp)
    swl1_ave(:)  = swl1_ave(:)  + swl1(:)     *(options%dtime/60._wp/60._wp/24._wp)
    swl2_ave(:)  = swl2_ave(:)  + swl2(:)     *(options%dtime/60._wp/60._wp/24._wp)
    swl3_ave(:)  = swl3_ave(:)  + swl3(:)     *(options%dtime/60._wp/60._wp/24._wp)
    swl4_ave(:)  = swl4_ave(:)  + swl4(:)     *(options%dtime/60._wp/60._wp/24._wp)
    lai_ave(:)   = lai_ave(:)   + lai(:)      *(options%dtime/60._wp/60._wp/24._wp)

    IF (is_newday(options%current_datetime, options%dtime)) THEN ! only once per day

      !!!! These will be removed after the FUEL issue 
      CALL calc_fuel_jsbach(       &
        & c_acid_ag1_ta(:),        & ! in
        & c_water_ag1_ta(:),       & ! in
        & c_ethanol_ag1_ta(:),     & ! in
        & c_nonsoluble_ag1_ta(:),  & ! in
        & c_acid_ag2_ta(:),        & ! in
        & c_water_ag2_ta(:),       & ! in
        & c_ethanol_ag2_ta(:),     & ! in
        & c_nonsoluble_ag2_ta(:),  & ! in
        & fract_fpc_max(:),        & ! in
        & fuel(:)                  & ! out
        & )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !if (SUM(fuel) > 0) then
      !  print*, SUM(fuel)
      !endif

      ! calculate auxiliary variables for DL fire model   <Rackhun Son>
      IF (dsl4jsb_Config(DISTURB_)%fire_algorithm == 5) THEN
        ! approximate daily climatology with 20-year window
        diff_year = iyr-time_config%tc_exp_startdate%date%year

        IF (is_time_experiment_start(options%current_datetime)) THEN   ! 1st iteration
          !t_air_clim(:,imo) = t_air_ave(:) /(options%dtime/60._wp/60._wp/24._wp)
          !q_clim(:,imo)     = q_ave(:)     /(options%dtime/60._wp/60._wp/24._wp)
          !q_rel_clim(:,imo) = q_rel_ave(:) /(options%dtime/60._wp/60._wp/24._wp)
          !!ws10_clim(:,imo)  = ws10_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          !swl1_clim(:,imo)  = swl1_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          !swl2_clim(:,imo)  = swl2_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          !swl3_clim(:,imo)  = swl3_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          !swl4_clim(:,imo)  = swl4_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          !lai_clim(:,imo)   = lai_ave(:)   /(options%dtime/60._wp/60._wp/24._wp)
          !fuel_yearly(:)    = fuel(:)

          t_air_clim(:,day_of_year) = t_air_ave(:) /(options%dtime/60._wp/60._wp/24._wp)
          q_clim(:,day_of_year)     = q_ave(:)     /(options%dtime/60._wp/60._wp/24._wp)
          q_rel_clim(:,day_of_year) = q_rel_ave(:) /(options%dtime/60._wp/60._wp/24._wp)
          swl1_clim(:,day_of_year)  = swl1_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          swl2_clim(:,day_of_year)  = swl2_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          swl3_clim(:,day_of_year)  = swl3_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          swl4_clim(:,day_of_year)  = swl4_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          lai_clim(:,day_of_year)   = lai_ave(:)   /(options%dtime/60._wp/60._wp/24._wp)
          fuel_yearly(:)            = fuel(:)

        !ELSEIF (iyr == time_config%tc_exp_startdate%date%year) THEN    ! 1st year
        ELSE

          !IF (iyr == time_config%tc_exp_startdate%date%year .AND. idy == 1) THEN   ! 1st day of month
          IF (iyr == time_config%tc_exp_startdate%date%year) THEN    ! 1st year
            !t_air_clim(:,imo) = t_air_ave(:)
            !q_clim(:,imo)     = q_ave(:)
            !q_rel_clim(:,imo) = q_rel_ave(:)
            !!ws10_clim(:,imo)  = ws10_ave(:)
            !swl1_clim(:,imo)  = swl1_ave(:)
            !swl2_clim(:,imo)  = swl2_ave(:)
            !swl3_clim(:,imo)  = swl3_ave(:)
            !swl4_clim(:,imo)  = swl4_ave(:)
            !lai_clim(:,imo)   = lai_ave(:)

            t_air_clim(:,day_of_year) = t_air_ave(:)
            q_clim(:,day_of_year)     = q_ave(:)
            q_rel_clim(:,day_of_year) = q_rel_ave(:)
            swl1_clim(:,day_of_year)  = swl1_ave(:)
            swl2_clim(:,day_of_year)  = swl2_ave(:)
            swl3_clim(:,day_of_year)  = swl3_ave(:)
            swl4_clim(:,day_of_year)  = swl4_ave(:)
            lai_clim(:,day_of_year)   = lai_ave(:)

          ELSEIF (day_of_year == 60) Then                                ! 2/29 <= 2/28
            t_air_clim(:,day_of_year) = t_air_clim(:,59)
            q_clim(:,day_of_year)     = q_clim(:,59)
            q_rel_clim(:,day_of_year) = q_rel_clim(:,59)
            swl1_clim(:,day_of_year)  = swl1_clim(:,59)
            swl2_clim(:,day_of_year)  = swl2_clim(:,59)
            swl3_clim(:,day_of_year)  = swl3_clim(:,59)
            swl4_clim(:,day_of_year)  = swl4_clim(:,59)
            lai_clim(:,day_of_year)   = lai_clim(:,59)

          ELSE
            !IF (diff_year < 1) Then
            IF (diff_year < 10) Then
              !clim_ratio = MAX(30._wp/idy, 1._wp)
              clim_ratio = MAX(10._wp/(diff_year+1), 1._wp)
            ELSE
              clim_ratio = 1._wp
            ENDIF
            !t_air_clim(:,imo) = t_air_clim(:,imo)*(30._wp-clim_ratio)/30._wp + t_air_ave(:)*clim_ratio/30._wp
            !q_clim(:,imo)     = q_clim(:,imo)*(30._wp-clim_ratio)/30._wp     + q_ave(:)*clim_ratio/30._wp
            !q_rel_clim(:,imo) = q_rel_clim(:,imo)*(30._wp-clim_ratio)/30._wp + q_rel_ave(:)*clim_ratio/30._wp
            !!ws10_clim(:,imo)  = ws10_clim(:,imo)*(30._wp-clim_ratio)/30._wp  + ws10_ave(:)*clim_ratio/30._wp
            !swl1_clim(:,imo)  = swl1_clim(:,imo)*(30._wp-clim_ratio)/30._wp  + swl1_ave(:)*clim_ratio/30._wp
            !swl2_clim(:,imo)  = swl2_clim(:,imo)*(30._wp-clim_ratio)/30._wp  + swl2_ave(:)*clim_ratio/30._wp
            !swl3_clim(:,imo)  = swl3_clim(:,imo)*(30._wp-clim_ratio)/30._wp  + swl3_ave(:)*clim_ratio/30._wp
            !swl4_clim(:,imo)  = swl4_clim(:,imo)*(30._wp-clim_ratio)/30._wp  + swl4_ave(:)*clim_ratio/30._wp
            !lai_clim(:,imo)   = lai_clim(:,imo)*(30._wp-clim_ratio)/30._wp   + lai_ave(:)*clim_ratio/30._wp

            t_air_clim(:,day_of_year) = t_air_clim(:,day_of_year)*(10._wp-clim_ratio)/10._wp + t_air_ave(:)*clim_ratio/10._wp
            q_clim(:,day_of_year)     = q_clim(:,day_of_year)*(10._wp-clim_ratio)/10._wp     + q_ave(:)*clim_ratio/10._wp
            q_rel_clim(:,day_of_year) = q_rel_clim(:,day_of_year)*(10._wp-clim_ratio)/10._wp + q_rel_ave(:)*clim_ratio/10._wp
            swl1_clim(:,day_of_year)  = swl1_clim(:,day_of_year)*(10._wp-clim_ratio)/10._wp  + swl1_ave(:)*clim_ratio/10._wp
            swl2_clim(:,day_of_year)  = swl2_clim(:,day_of_year)*(10._wp-clim_ratio)/10._wp  + swl2_ave(:)*clim_ratio/10._wp
            swl3_clim(:,day_of_year)  = swl3_clim(:,day_of_year)*(10._wp-clim_ratio)/10._wp  + swl3_ave(:)*clim_ratio/10._wp
            swl4_clim(:,day_of_year)  = swl4_clim(:,day_of_year)*(10._wp-clim_ratio)/10._wp  + swl4_ave(:)*clim_ratio/10._wp
            lai_clim(:,day_of_year)   = lai_clim(:,day_of_year)*(10._wp-clim_ratio)/10._wp   + lai_ave(:)*clim_ratio/10._wp
          ENDIF

          IF (iyr == time_config%tc_exp_startdate%date%year) THEN   ! 1st year
            yearly_ratio = MAX(365._wp/get_year_day(options%current_datetime), 1._wp)
          ELSE
            yearly_ratio = 1._wp
          ENDIF
          fuel_yearly(:)= fuel_yearly(:)*(365._wp-yearly_ratio)/365._wp  + fuel(:)*yearly_ratio/365._wp

        ENDIF

        !if (isnan(fuel(1))) then
        !  print*, 'fuel'
        !endif

        ! modify pft fractions for Fire
        !modified_pfts_fract(:,1) = fract_snow_soil(:)
        !modified_pfts_fract(:,2) = pfts_fract(:,1)*(1-fract_snow_soil(:))*fract_fpc_max     ! Tropical broadleaf evergreen
        !modified_pfts_fract(:,3) = pfts_fract(:,2)*(1-fract_snow_soil(:))*fract_fpc_max     ! Tropical broadleaf deciduous
        !modified_pfts_fract(:,4) = pfts_fract(:,3)*(1-fract_snow_soil(:))*fract_fpc_max     ! Extra-tropical broadleaf evergreen
        !modified_pfts_fract(:,5) = pfts_fract(:,4)*(1-fract_snow_soil(:))*fract_fpc_max     ! Extra_tropical broadleaf deciduous
        !modified_pfts_fract(:,6) = pfts_fract(:,5)*(1-fract_snow_soil(:))*fract_fpc_max     ! Raingreen Shrubs
        !modified_pfts_fract(:,7) = pfts_fract(:,6)*(1-fract_snow_soil(:))*fract_fpc_max     ! Deciduous Shrubs
        !modified_pfts_fract(:,8) = SUM(pfts_fract(:,7:8),DIM=2)*(1-fract_snow_soil(:))*fract_fpc_max   ! Grass
        !modified_pfts_fract(:,9) = (1-fract_snow_soil(:))*(1-fract_fpc_max)                 ! Others

        modified_pfts_fract(:,1) = fract_snow(:)
        modified_pfts_fract(:,2) = pfts_fract(:,1)*(1-fract_snow(:))*fract_fpc_max     ! Tropical broadleaf evergreen
        modified_pfts_fract(:,3) = pfts_fract(:,2)*(1-fract_snow(:))*fract_fpc_max     ! Tropical broadleaf deciduous
        modified_pfts_fract(:,4) = pfts_fract(:,3)*(1-fract_snow(:))*fract_fpc_max     ! Extra-tropical broadleaf evergreen
        modified_pfts_fract(:,5) = pfts_fract(:,4)*(1-fract_snow(:))*fract_fpc_max     ! Extra_tropical broadleaf deciduous
        modified_pfts_fract(:,6) = pfts_fract(:,5)*(1-fract_snow(:))*fract_fpc_max     ! Raingreen Shrubs
        modified_pfts_fract(:,7) = pfts_fract(:,6)*(1-fract_snow(:))*fract_fpc_max     ! Deciduous Shrubs
        modified_pfts_fract(:,8) = SUM(pfts_fract(:,7:8),DIM=2)*(1-fract_snow(:))*fract_fpc_max   ! Grass
        modified_pfts_fract(:,9) = (1-fract_snow(:))*(1-fract_fpc_max)                 ! Others

        lightning_climatology => dsl4jsb_memory(DISTURB_)%lightning_climatology(ics:ice, iblk, day_of_year)    ! <Rackhun Son>
        population_density    => dsl4jsb_memory(DISTURB_)%population_density(ics:ice, iblk, day_of_year)       ! <Rackhun Son>

        t_air_local_clim      => dsl4jsb_memory(DISTURB_)%t_air_local_clim(ics:ice, iblk, day_of_year)         ! <Rackhun Son>
        q_rel_local_clim      => dsl4jsb_memory(DISTURB_)%q_rel_local_clim(ics:ice, iblk, day_of_year)         ! <Rackhun Son>
        q_local_clim          => dsl4jsb_memory(DISTURB_)%q_local_clim(ics:ice, iblk, day_of_year)             ! <Rackhun Son>
        ws10_local_clim       => dsl4jsb_memory(DISTURB_)%ws10_local_clim(ics:ice, iblk, day_of_year)          ! <Rackhun Son>
        swvl1_local_clim      => dsl4jsb_memory(DISTURB_)%swvl1_local_clim(ics:ice, iblk, day_of_year)         ! <Rackhun Son>
        swvl2_local_clim      => dsl4jsb_memory(DISTURB_)%swvl2_local_clim(ics:ice, iblk, day_of_year)         ! <Rackhun Son>
        swvl3_local_clim      => dsl4jsb_memory(DISTURB_)%swvl3_local_clim(ics:ice, iblk, day_of_year)         ! <Rackhun Son>
        swvl4_local_clim      => dsl4jsb_memory(DISTURB_)%swvl4_local_clim(ics:ice, iblk, day_of_year)         ! <Rackhun Son>
        lai_local_clim        => dsl4jsb_memory(DISTURB_)%lai_local_clim(ics:ice, iblk, day_of_year)           ! <Rackhun Son>
        !t_air_local_clim      => dsl4jsb_memory(DISTURB_)%t_air_local_clim(ics:ice, iblk, imo)         ! <Rackhun Son>
        !q_rel_local_clim      => dsl4jsb_memory(DISTURB_)%q_rel_local_clim(ics:ice, iblk, imo)         ! <Rackhun Son>
        !q_local_clim          => dsl4jsb_memory(DISTURB_)%q_local_clim(ics:ice, iblk, imo)             ! <Rackhun Son>
        !!ws10_local_clim       => dsl4jsb_memory(DISTURB_)%ws10_local_clim(ics:ice, iblk, imo)          ! <Rackhun Son>
        !swvl1_local_clim      => dsl4jsb_memory(DISTURB_)%swvl1_local_clim(ics:ice, iblk, imo)         ! <Rackhun Son>
        !swvl2_local_clim      => dsl4jsb_memory(DISTURB_)%swvl2_local_clim(ics:ice, iblk, imo)         ! <Rackhun Son>
        !swvl3_local_clim      => dsl4jsb_memory(DISTURB_)%swvl3_local_clim(ics:ice, iblk, imo)         ! <Rackhun Son>
        !swvl4_local_clim      => dsl4jsb_memory(DISTURB_)%swvl4_local_clim(ics:ice, iblk, imo)         ! <Rackhun Son>
        !lai_local_clim        => dsl4jsb_memory(DISTURB_)%lai_local_clim(ics:ice, iblk, imo)           ! <Rackhun Son>
        
      ENDIF

      !print*, t_air_min_bias(1)

      ! reset daily statistics
      t_air_ave(:) = 0._wp
      q_ave(:) = 0._wp
      q_rel_ave(:) = 0._wp
      swl1_ave(:) = 0._wp
      swl2_ave(:) = 0._wp
      swl3_ave(:) = 0._wp
      swl4_ave(:) = 0._wp
      lai_ave(:) = 0._wp

      burned_fract(:) = 0._wp
      IF (      dsl4jsb_Lctlib_param(dynamic_PFT) &
          .OR. (dsl4jsb_Config(DISTURB_)%lburn_pasture .AND. dsl4jsb_Lctlib_param(pasture_PFT)) &
         ) THEN

        SELECT CASE (dsl4jsb_Config(DISTURB_)%fire_algorithm)
          CASE (0) !! No fire algorithm
          CASE (1) !! jsbach algorithm

            ! Calculate burned area. Woody and grass types call the subroutine with different parameters
            IF (dsl4jsb_Lctlib_param(woody_PFT)) THEN  ! Woody Type
              fire_minimum = dsl4jsb_Config(DISTURB_)%fire_minimum_woody
              fire_tau     = dsl4jsb_Config(DISTURB_)%fire_tau_woody
            ELSE                             ! all other types
              fire_minimum = dsl4jsb_Config(DISTURB_)%fire_minimum_grass
              fire_tau     = dsl4jsb_Config(DISTURB_)%fire_tau_grass
            ENDIF

            CALL burned_fract_jsbach(                            &
              & dsl4jsb_Config(DISTURB_)%fire_rel_hum_threshold, &
              & dsl4jsb_Config(DISTURB_)%fire_litter_threshold,  &
              & fire_minimum,                                    & ! in
              & fire_tau,                                        & ! in
              & q_rel_air_climbuf(:),                            & ! in
              & fuel(:),                                         & ! in
              & burned_fract(:)                                  & ! inout
              & )
          CASE (2) !! Arora & Boer algorithm
            CALL finish('disturbed_fract','Arora & Boer algorithm not implemented yet.')
          CASE (3)
            CALL finish('disturbed_fract','Thornicke spit fire not implemented yet.')
          CASE (4) !! Read burned_fract
          CASE (5) !! DL based new fire algoritm   <Rackhun Son>
            !IF (maxval(fuel) /= 0.) print*, minval(fuel), maxval(fuel)
            !if (SUM(c_acid_ag1_ta)>0) then
            !  print*, SUM(c_acid_ag1_ta)
            !endif

            CALL burned_fract_dl(                                      &
              & t_air(:),                                              &
              & t_air_local_clim(:), t_air_clim(:,day_of_year),        &  
              & q_rel_air(:),                                          &
              & q_rel_local_clim(:), q_rel_clim(:,day_of_year),        &
              & q_air(:),                                              &
              & q_local_clim(:),     q_clim(:,day_of_year),            &
              & wind_10m(:),                                           &
              & rain(:),                                               &
              & swl1(:),                                               &
              & swvl1_local_clim(:), swl1_clim(:,day_of_year),         &
              & swl2(:),                                               &
              & swvl2_local_clim(:), swl2_clim(:,day_of_year),         &
              & swl3(:),                                               &
              & swvl3_local_clim(:), swl3_clim(:,day_of_year),         &
              & swl4(:),                                               &
              & swvl4_local_clim(:), swl4_clim(:,day_of_year),         &
              & lai(:),                                                &
              & lai_local_clim(:),   lai_clim(:,day_of_year),          &
              & modified_pfts_fract(:,:),                              &
              & fuel(:),             fuel_yearly(:),                   &
             !& fuel(:),                                               &
              & population_density(:),                                 &
              & lightning_climatology(:),                              &
              & elevation(:),                                          &

              & nc,                                                    &
              & tile%name,                                             &

              & dsl4jsb_Config(DISTURB_)%d_climate_input,              &
              & dsl4jsb_Config(DISTURB_)%d_climate_hidden,             &
              & dsl4jsb_Config(DISTURB_)%d_vegetation_input,           &
              & dsl4jsb_Config(DISTURB_)%d_vegetation_hidden,          &
              & dsl4jsb_Config(DISTURB_)%d_static_input,               &
              & dsl4jsb_Config(DISTURB_)%d_static_hidden,              &
              & dsl4jsb_Config(DISTURB_)%d_output,                     &
              & dsl4jsb_Config(DISTURB_)%d_emb_input,                  &
              & dsl4jsb_Config(DISTURB_)%d_emb_hidden,                 &
              & dsl4jsb_Config(DISTURB_)%d_emb_output,                 &

              & dsl4jsb_Config(DISTURB_)%climate_i2i_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_i2i_bias,             &
              & dsl4jsb_Config(DISTURB_)%climate_h2i_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_h2i_bias,             &
              & dsl4jsb_Config(DISTURB_)%climate_i2f_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_i2f_bias,             &
              & dsl4jsb_Config(DISTURB_)%climate_h2f_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_h2f_bias,             &
              & dsl4jsb_Config(DISTURB_)%climate_i2g_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_i2g_bias,             &
              & dsl4jsb_Config(DISTURB_)%climate_h2g_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_h2g_bias,             &
              & dsl4jsb_Config(DISTURB_)%climate_i2o_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_i2o_bias,             &
              & dsl4jsb_Config(DISTURB_)%climate_h2o_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_h2o_bias,             &
              & dsl4jsb_Config(DISTURB_)%climate_out_weight,           &
              & dsl4jsb_Config(DISTURB_)%climate_out_bias,             &

              & dsl4jsb_Config(DISTURB_)%vegetation_i2i_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_i2i_bias,          &
              & dsl4jsb_Config(DISTURB_)%vegetation_h2i_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_h2i_bias,          &
              & dsl4jsb_Config(DISTURB_)%vegetation_i2f_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_i2f_bias,          &
              & dsl4jsb_Config(DISTURB_)%vegetation_h2f_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_h2f_bias,          &
              & dsl4jsb_Config(DISTURB_)%vegetation_i2g_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_i2g_bias,          &
              & dsl4jsb_Config(DISTURB_)%vegetation_h2g_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_h2g_bias,          &
              & dsl4jsb_Config(DISTURB_)%vegetation_i2o_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_i2o_bias,          &
              & dsl4jsb_Config(DISTURB_)%vegetation_h2o_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_h2o_bias,          &
              & dsl4jsb_Config(DISTURB_)%vegetation_out_weight,        &
              & dsl4jsb_Config(DISTURB_)%vegetation_out_bias,          &

              & dsl4jsb_Config(DISTURB_)%static_layer1_weight,         &
              & dsl4jsb_Config(DISTURB_)%static_layer1_bias,           &
              & dsl4jsb_Config(DISTURB_)%static_layer2_weight,         &
              & dsl4jsb_Config(DISTURB_)%static_layer2_bias,           &

              & dsl4jsb_Config(DISTURB_)%climate_emb_layer1_weight,    &
              & dsl4jsb_Config(DISTURB_)%climate_emb_layer1_bias,      &
              & dsl4jsb_Config(DISTURB_)%climate_emb_layer2_weight,    &
              & dsl4jsb_Config(DISTURB_)%climate_emb_layer2_bias,      &

              & dsl4jsb_Config(DISTURB_)%vegetation_emb_layer1_weight, &
              & dsl4jsb_Config(DISTURB_)%vegetation_emb_layer1_bias,   &
              & dsl4jsb_Config(DISTURB_)%vegetation_emb_layer2_weight, &
              & dsl4jsb_Config(DISTURB_)%vegetation_emb_layer2_bias,   &

              & dsl4jsb_Config(DISTURB_)%t_air_minmax,                 &
              & dsl4jsb_Config(DISTURB_)%t_air_ano_minmax,             &
              & dsl4jsb_Config(DISTURB_)%t_air_clim_minmax,            &
              & dsl4jsb_Config(DISTURB_)%q_rel_minmax,                 &
              & dsl4jsb_Config(DISTURB_)%q_rel_ano_minmax,             &
              & dsl4jsb_Config(DISTURB_)%q_rel_clim_minmax,            &
              & dsl4jsb_Config(DISTURB_)%q_minmax,                     &
              & dsl4jsb_Config(DISTURB_)%q_ano_minmax,                 &
              & dsl4jsb_Config(DISTURB_)%q_clim_minmax,                &
              & dsl4jsb_Config(DISTURB_)%wind_10m_minmax,              &
              & dsl4jsb_Config(DISTURB_)%rain_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%swl1_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%swl1_ano_minmax,              &
              & dsl4jsb_Config(DISTURB_)%swl1_clim_minmax,             &
              & dsl4jsb_Config(DISTURB_)%swl2_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%swl2_ano_minmax,              &
              & dsl4jsb_Config(DISTURB_)%swl2_clim_minmax,             &
              & dsl4jsb_Config(DISTURB_)%swl3_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%swl3_ano_minmax,              &
              & dsl4jsb_Config(DISTURB_)%swl3_clim_minmax,             &
              & dsl4jsb_Config(DISTURB_)%swl4_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%swl4_ano_minmax,              &
              & dsl4jsb_Config(DISTURB_)%swl4_clim_minmax,             &
              & dsl4jsb_Config(DISTURB_)%lai_minmax,                   &
              & dsl4jsb_Config(DISTURB_)%lai_ano_minmax,               &
              & dsl4jsb_Config(DISTURB_)%lai_clim_minmax,              &
              & dsl4jsb_Config(DISTURB_)%population_density_minmax,    &
              & dsl4jsb_Config(DISTURB_)%lightning_climatology_minmax, &
              & dsl4jsb_Config(DISTURB_)%elevation_minmax,             &
              & dsl4jsb_Config(DISTURB_)%agb_minmax,                   &
              & dsl4jsb_Config(DISTURB_)%pft1_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%pft2_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%pft3_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%pft4_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%pft5_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%pft6_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%pft7_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%pft8_minmax,                  &
              & dsl4jsb_Config(DISTURB_)%pft9_minmax,                  &
              
              & climate_lstm_hidden(:,:),                              & ! inout
              & climate_lstm_cell(:,:),                                & ! inout
              & vegetation_lstm_hidden(:,:),                           & ! inout
              & vegetation_lstm_cell(:,:),                             & ! inout
              & burned_fract(:)                                        & ! inout
              & )

          CASE DEFAULT
            CALL finish('disturbed_fract','Unknown fire algorithm')
        END SELECT
      END IF

      damaged_fract(:) = 0._wp
      IF (dsl4jsb_Lctlib_param(woody_PFT) .AND. dsl4jsb_Lctlib_param(dynamic_PFT)) THEN
        SELECT CASE (dsl4jsb_Config(DISTURB_)%windbreak_algorithm)
          CASE (0) !! No windbreak algorithm
          CASE (1) !! jsbach algorithm
            CALL broken_woody_fract_jsbach(                 & ! in
              & dsl4jsb_Config(DISTURB_)%wnd_threshold,     & ! in
              & dsl4jsb_Config(DISTURB_)%wnd_damage_scale,  & ! in
              & cover_fract_pot(:),                         & ! in
              & prev_day_max_wind_10m(:),                   & ! in
              & max_wind_10m(:),                            & ! in
              & damaged_fract(:)                            & ! inout
              & )
          CASE DEFAULT
            CALL finish('disturbed_frac','Unknown windbreak algorithm')
        END SELECT
      ENDIF

      ! Init diagnostic carbon fluxes
      IF ( tile%Has_process_memory(CARBON_) ) THEN
        co2flux_fire_all_2_atm(:) = 0._wp
        cflux_dist_green_2_soil(:) = 0._wp
        cflux_dist_woods_2_soil(:) = 0._wp
      ENDIF

    ENDIF ! is_newday

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_damaged_area
  
  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "damaged_area"
  !
  SUBROUTINE aggregate_damaged_area(tile, options)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DISTURB_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_damaged_area'

    INTEGER  :: iblk, ics, ice
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(DISTURB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(DISTURB_, burned_fract,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, damaged_fract,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, q_rel_air_climbuf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, q_rel_air_climbuf_yDay,  weighted_by_fract)

  END SUBROUTINE aggregate_damaged_area

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "natural_disturbances"
  !>
  !> Task "natural_disturbances" calculates the disturbance of plants due to fire and windbreak
  !> if flcc and wlcc are not active, damaged area and carbon effects are calculated on this tile
  !> else only ta variables are updated
  !
  SUBROUTINE update_natural_disturbances(tile, options)

    USE mo_disturb_process,       ONLY: burned_fract_jsbach, broken_woody_fract_jsbach
    USE mo_carbon_process,        ONLY: relocate_carbon_fire, relocate_carbon_damage
    USE mo_carbon_interface,      ONLY: recalc_per_tile_vars, calculate_current_c_ag_1_and_bg_sums,  &
                                      & calculate_current_c_ta_state_sum, check_carbon_conservation, &
                                      & recalc_carbon_per_tile_vars
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_time,              ONLY: is_time_experiment_start, timesteps_per_day, is_newday

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    REAL(wp), DIMENSION(options%nc)           :: old_c_state_sum_ta, dummy_flux

    INTEGER :: iblk, nc, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_natural_disturbances'

    ! Declare process memories
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_config(DISTURB_)
    dsl4jsb_Def_config(FLCC_)
    dsl4jsb_Def_config(WLCC_)
    dsl4jsb_Def_config(CARBON_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: cconservation_fire
    dsl4jsb_Real2D_onChunk :: cconservation_wind
    !
    dsl4jsb_Real2D_onChunk :: burned_fract
    dsl4jsb_Real2D_onChunk :: damaged_fract
    !
    dsl4jsb_Real2D_onChunk :: c_green
    dsl4jsb_Real2D_onChunk :: c_woods
    dsl4jsb_Real2D_onChunk :: c_reserve
    dsl4jsb_Real2D_onChunk :: cflux_dist_green_2_soil
    dsl4jsb_Real2D_onChunk :: cflux_dist_woods_2_soil
    !
    dsl4jsb_Real2D_onChunk :: c_acid_ag1
    dsl4jsb_Real2D_onChunk :: c_water_ag1
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag1
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag1
    dsl4jsb_Real2D_onChunk :: c_acid_ag2
    dsl4jsb_Real2D_onChunk :: c_water_ag2
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag2
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag2
    ! below ground
    dsl4jsb_Real2D_onChunk :: c_acid_bg1
    dsl4jsb_Real2D_onChunk :: c_water_bg1
    dsl4jsb_Real2D_onChunk :: c_ethanol_bg1
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_bg1
    dsl4jsb_Real2D_onChunk :: c_acid_bg2
    dsl4jsb_Real2D_onChunk :: c_water_bg2
    dsl4jsb_Real2D_onChunk :: c_ethanol_bg2
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_bg2
    !
    dsl4jsb_Real2D_onChunk :: c_humus_1
    dsl4jsb_Real2D_onChunk :: c_humus_2
    dsl4jsb_Real2D_onChunk :: co2flux_fire_all_2_atm
    !
    dsl4jsb_Real2D_onChunk :: co2flux_fire_all_2_atm_ta
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(DISTURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process configurations
    dsl4jsb_Get_config(DISTURB_)
    dsl4jsb_Get_config(CARBON_)

    ! Get process memories
    dsl4jsb_Get_memory(DISTURB_)
    dsl4jsb_Get_memory(CARBON_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   burned_fract)             ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   damaged_fract)            ! out
    !
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   cconservation_wind)     ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   cconservation_fire)     ! out
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_green)                  ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_woods)                  ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_reserve)                ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_green_2_soil)  ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_woods_2_soil)  ! inout
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_ag1)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_ag1)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_ag1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_ag1)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_ag2)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_ag2)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_ag2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_ag2)        ! inout
    ! below ground
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_bg1)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_bg1)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_bg1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_bg1)        ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_bg2)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_bg2)             ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_bg2)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_bg2)        ! inout
    !
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_humus_1)               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_humus_2)               ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    co2flux_fire_all_2_atm)  ! inout

    !ta variable required for the carbon conservation test
    dsl4jsb_Get_var2D_onChunk(CARBON_,    co2flux_fire_all_2_atm_ta) ! out

    ! -------------------------------------------------------------------------------------------------- !
    IF (is_newday(options%current_datetime, options%dtime)) THEN ! only once per day

      IF ( tile%Has_process_memory(CARBON_) ) THEN

        ! Calculate carbon relocation connected with damages and fires
        IF (.NOT. model%processes(WLCC_)%p%config%active) THEN
          ! If WLCC is not active, wind disturbance is executed directly on this pft
          ! (If WLCC would have been active, the carbon consequences of wind disturbance would have already been calculated)
          SELECT CASE (dsl4jsb_Config(DISTURB_)%windbreak_algorithm)
            CASE (0) !! No windbreak algorithm
            CASE (1) !! jsbach algorithm

              !Calculate current sum before operation for c conservation test
              CALL recalc_carbon_per_tile_vars(tile, options)
              CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum_ta(:))

              CALL relocate_carbon_damage(                 &
                & damaged_fract(:),                        & ! in
                & dsl4jsb_Lctlib_param(LeafLit_coef(1:5)), & ! in
                & dsl4jsb_Lctlib_param(WoodLit_coef(1:5)), & ! in
                & c_green(:),                              & ! inout
                & c_woods(:),                              & ! inout
                & c_reserve(:),                            & ! inout
                & cflux_dist_green_2_soil(:),              & ! inout
                & cflux_dist_woods_2_soil(:),              & ! inout
                & c_acid_ag1(:),                           & ! inout
                & c_water_ag1(:),                          & ! inout
                & c_ethanol_ag1(:),                        & ! inout
                & c_nonsoluble_ag1(:),                     & ! inout
                & c_acid_ag2(:),                           & ! inout
                & c_water_ag2(:),                          & ! inout
                & c_ethanol_ag2(:),                        & ! inout
                & c_nonsoluble_ag2(:),                     & ! inout
                & c_acid_bg1(:),                           & ! inout
                & c_water_bg1(:),                          & ! inout
                & c_ethanol_bg1(:),                        & ! inout
                & c_nonsoluble_bg1(:),                     & ! inout
                & c_acid_bg2(:),                           & ! inout
                & c_water_bg2(:),                          & ! inout
                & c_ethanol_bg2(:),                        & ! inout
                & c_nonsoluble_bg2(:),                     & ! inout
                & c_humus_1(:),                            & ! inout
                & c_humus_2(:)                             & ! inout
                & )

              CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
              CALL recalc_per_tile_vars(tile, options,                                      &
                & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
                & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
                & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
                & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
                & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
                & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
                & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
                & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
                & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
                & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
                & cflux_dist_green_2_soil = cflux_dist_green_2_soil(:),                     &
                & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:))

              ! wind has no extra flux to atmos, only redistribution among cpools - therefore cflux = 0._wp
              dummy_flux = 0._wp
              CALL check_carbon_conservation(tile, options, old_c_state_sum_ta(:), &
                & dummy_flux(:), cconservation_wind(:))

            CASE DEFAULT
              CALL finish('disturbed_fract','Unknown windbreak algorithm')
          END SELECT
        ENDIF

        IF (.NOT. model%processes(FLCC_)%p%config%active) THEN
          ! If FLCC is not active, fire disturbance is executed directly on this pft
          ! (If FLCC would have been active, the carbon consequences of fire disturbance would have already been calculated)

          co2flux_fire_all_2_atm = 0._wp
          SELECT CASE (dsl4jsb_Config(DISTURB_)%fire_algorithm)
            CASE (0) !! No fire algorithm
            CASE (1) !! jsbach algorithm

              !Calculate current sum before operation for c conservation test
              CALL recalc_carbon_per_tile_vars(tile, options)
              CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum_ta(:))

              CALL relocate_carbon_fire(                    &
                & burned_fract(:),                          & ! in, only for windbreak
                & dsl4jsb_Config(CARBON_)%fire_fract_wood_2_atmos , & ! in
                & dsl4jsb_Lctlib_param(LeafLit_coef(1:5)),  & ! in
                & dsl4jsb_Lctlib_param(WoodLit_coef(1:5)),  & ! in
                & c_green(:),                               & ! inout
                & c_woods(:),                               & ! inout
                & c_reserve(:),                             & ! inout
                & cflux_dist_green_2_soil(:),               & ! inout
                & cflux_dist_woods_2_soil(:),               & ! inout
                & c_acid_ag1(:),                            & ! inout
                & c_water_ag1(:),                           & ! inout
                & c_ethanol_ag1(:),                         & ! inout
                & c_nonsoluble_ag1(:),                      & ! inout
                & c_acid_ag2(:),                            & ! inout
                & c_water_ag2(:),                           & ! inout
                & c_ethanol_ag2(:),                         & ! inout
                & c_nonsoluble_ag2(:),                      & ! inout
                & c_acid_bg1(:),                            & ! inout
                & c_water_bg1(:),                           & ! inout
                & c_ethanol_bg1(:),                         & ! inout
                & c_nonsoluble_bg1(:),                      & ! inout
                & c_acid_bg2(:),                            & ! inout
                & c_water_bg2(:),                           & ! inout
                & c_ethanol_bg2(:),                         & ! inout
                & c_nonsoluble_bg2(:),                      & ! inout
                & c_humus_1(:),                             & ! inout
                & c_humus_2(:),                             & ! inout
                & co2flux_fire_all_2_atm                    & ! out
                & )

              CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
              CALL recalc_per_tile_vars(tile, options,                                      &
                & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
                & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
                & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
                & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
                & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
                & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
                & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
                & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
                & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
                & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
                & cflux_dist_green_2_soil = cflux_dist_green_2_soil(:),                     &
                & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:),                     &
                & co2flux_fire_all_2_atm = co2flux_fire_all_2_atm(:))

              ! For conservation test: negate currently positve co2flux_fire_all_2_atm_ta (flux away from land!)
              ! and convert from CO2 flux to carbon change per day

              CALL check_carbon_conservation(tile, options, old_c_state_sum_ta(:), &
                & -co2flux_fire_all_2_atm_ta(:) * sec_per_day / molarMassCO2_kg, cconservation_fire(:))

            CASE (2) !! Arora & Boer algorithm
              CALL finish('disturbed_fract','Arora & Boer algorithm not implemented yet.')
            CASE (3)
              CALL finish('disturbed_fract','Thornicke spit fire not implemented yet.')
            CASE (4) !! Read burned_fract
            CASE (5) !! jsbach algorithm

              !Calculate current sum before operation for c conservation test
              CALL recalc_carbon_per_tile_vars(tile, options)
              CALL calculate_current_c_ta_state_sum(tile, options, old_c_state_sum_ta(:))

              CALL relocate_carbon_fire(                    &
                & burned_fract(:),                          & ! in, only for windbreak
                & dsl4jsb_Config(CARBON_)%fire_fract_wood_2_atmos , & ! in
                & dsl4jsb_Lctlib_param(LeafLit_coef(1:5)),  & ! in
                & dsl4jsb_Lctlib_param(WoodLit_coef(1:5)),  & ! in
                & c_green(:),                               & ! inout
                & c_woods(:),                               & ! inout
                & c_reserve(:),                             & ! inout
                & cflux_dist_green_2_soil(:),               & ! inout
                & cflux_dist_woods_2_soil(:),               & ! inout
                & c_acid_ag1(:),                            & ! inout
                & c_water_ag1(:),                           & ! inout
                & c_ethanol_ag1(:),                         & ! inout
                & c_nonsoluble_ag1(:),                      & ! inout
                & c_acid_ag2(:),                            & ! inout
                & c_water_ag2(:),                           & ! inout
                & c_ethanol_ag2(:),                         & ! inout
                & c_nonsoluble_ag2(:),                      & ! inout
                & c_acid_bg1(:),                            & ! inout
                & c_water_bg1(:),                           & ! inout
                & c_ethanol_bg1(:),                         & ! inout
                & c_nonsoluble_bg1(:),                      & ! inout
                & c_acid_bg2(:),                            & ! inout
                & c_water_bg2(:),                           & ! inout
                & c_ethanol_bg2(:),                         & ! inout
                & c_nonsoluble_bg2(:),                      & ! inout
                & c_humus_1(:),                             & ! inout
                & c_humus_2(:),                             & ! inout
                & co2flux_fire_all_2_atm                    & ! out
                & )

              CALL calculate_current_c_ag_1_and_bg_sums(tile, options)
              CALL recalc_per_tile_vars(tile, options,                                      &
                & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:),      &
                & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:),                 &
                & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
                & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:),                 &
                & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
                & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:),                 &
                & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
                & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:),                 &
                & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
                & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:),                       &
                & cflux_dist_green_2_soil = cflux_dist_green_2_soil(:),                     &
                & cflux_dist_woods_2_soil = cflux_dist_woods_2_soil(:),                     &
                & co2flux_fire_all_2_atm = co2flux_fire_all_2_atm(:))

              ! For conservation test: negate currently positve co2flux_fire_all_2_atm_ta (flux away from land!)
              ! and convert from CO2 flux to carbon change per day

              CALL check_carbon_conservation(tile, options, old_c_state_sum_ta(:), &
                & -co2flux_fire_all_2_atm_ta(:) * sec_per_day / molarMassCO2_kg, cconservation_fire(:))
            CASE DEFAULT
              CALL finish('disturbed_fract','Unknown fire algorithm')
          END SELECT
        ENDIF
      ENDIF ! CARBON active

    ENDIF ! is_newday

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_natural_disturbances

  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "natural_disturbances"
  !
  ! -------------------------------------------------------------------------------------------------------!
  SUBROUTINE aggregate_natural_disturbances(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(CARBON_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_natural_disturbance'

    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(CARBON_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    IF ( tile%Is_process_calculated(CARBON_) ) THEN
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_green_ta,              weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_woods_ta,              weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_reserve_ta,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_green_2_soil_ta, weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta, weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag1_ta,           weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag1_ta,          weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag1_ta,        weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag1_ta,     weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag2_ta,           weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag2_ta,          weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag2_ta,        weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag2_ta,     weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_1_ta,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_2_ta,            weighted_by_fract)
       dsl4jsb_Aggregate_onChunk(CARBON_,  co2flux_fire_all_2_atm_ta, weighted_by_fract)
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_natural_disturbances

#endif
END MODULE mo_disturb_interface
