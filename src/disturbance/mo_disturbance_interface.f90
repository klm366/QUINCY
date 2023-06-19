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
  dsl4jsb_Use_processes PHENO_, HYDRO_, ANTHRO2L_   ! <Rackhun Son>

  ! Use of process configurations (t_disturb_config)
  dsl4jsb_Use_config(DISTURB_)
  dsl4jsb_Use_config(CARBON_)

  ! Use of process memories (t_disturb_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(DISTURB_)
  dsl4jsb_Use_memory(FUEL_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(PHENO_)       ! <Rackhun Son>
  dsl4jsb_Use_memory(HYDRO_)       ! <Rackhun Son>
  dsl4jsb_Use_memory(ANTHRO2L_)    ! <Rackhun Son>

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
    USE mo_disturb_process,       ONLY: burned_fract_dl, get_relative_humidity_air_v2, Koppen_Geiger_classification_lv2 ! <Rackhun Son>  
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_time,              ONLY: is_time_experiment_start, timesteps_per_day, is_newday
    USE mo_disturbance_constants, ONLY: persist_rel_hum, persist_wind_10m

    USE mo_jsb_time_iface,        ONLY: get_year_day    ! <Rackhun Son>
    USE mo_time_config,           ONLY: time_config     ! <Rackhun Son>
    USE mo_jsb_grid,              ONLY: Get_grid        ! <Rackhun Son>
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid      ! <Rackhun Son>
    USE mo_global_variable                              ! <Rackhun Son>

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    TYPE(t_jsb_grid),  POINTER                :: grid  ! <Rackhun Son>
    INTEGER :: iblk, nc, ics, ice
    INTEGER :: iyr, imo, idy, day_of_year, i           ! <Rackhun Son>
    REAL(wp) :: clim_ratio, yearly_ratio, diff_year    ! average ratio for climatology <Rackhun Son>

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
    dsl4jsb_Def_memory(ANTHRO2L_)                    ! <Rackhun Son>

    dsl4jsb_Def_config(DISTURB_)

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

    ! axuxiliary input variables for fire
    dsl4jsb_Real2D_onChunk :: t_air                  ! temperature [K]   <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: rain                   ! precipitation <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: lai                    ! LAI <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: fract_fpc_max          ! Maximum fraction of tile that can be covered by leave <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: pfts_fract             ! fraction of PFTs <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: fract_snow             ! Snow fraction of tile <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: elevation              ! elevation <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl            ! water content in soil layers [m] <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: soil_depth_sl          ! soil depth within each soil layer [m] <Rackhun Son>

    dsl4jsb_Real2D_onChunk :: t_air_ave              ! daily average of temperature <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: q_ave                  ! daily average of specific humidity <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: q_rel_ave              ! daily average of relative humidity <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: ws10_ave               ! daily average of wind speed <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: rain_ave               ! daily average of precipiataion <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl1_ave               ! daily average of volumetric soil water content in the 1st layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl2_ave               ! daily average of volumetric soil water content in the 2nd layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl3_ave               ! daily average of volumetric soil water content in the 3rd layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl4_ave               ! daily average of volumetric soil water content in the 4th layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: swl5_ave               ! daily average of volumetric soil water content in the 4th layer <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: lai_ave                ! daily average of LAI <Rackhun Son>

    dsl4jsb_Real3D_onChunk :: t_air_annual           ! annual record of temperature <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: rain_annual            ! annual record of precipitation <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: t_air_clim             ! daily climatology of temperature <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: q_clim                 ! daily climatology of specific humidity <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: q_rel_clim             ! daily climatology of relative humidity <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl1_clim              ! daily climatology of volumetric soil water content in the 1st layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl2_clim              ! daily climatology of volumetric soil water content in the 2nd layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl3_clim              ! daily climatology of volumetric soil water content in the 3rd layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl4_clim              ! daily climatology of volumetric soil water content in the 4th layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: swl5_clim              ! daily climatology of volumetric soil water content in the 4th layer <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: lai_clim               ! daily climatology of LAI <Rackhun Son>

    dsl4jsb_Real3D_onChunk :: V1_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V2_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V3_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V4_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V5_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V6_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V7_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V8_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V9_memory              ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V10_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V11_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V12_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V13_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V14_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V15_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V16_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V17_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V18_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V19_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V20_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V21_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V22_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V23_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V24_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V25_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V26_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V27_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V28_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V29_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V30_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V31_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V32_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V33_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V34_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V35_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V36_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V37_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V38_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V39_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V40_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V41_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V42_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V43_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V44_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V45_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V46_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V47_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V48_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V49_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V50_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V51_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V52_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: V53_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: RHC_memory             ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: Y_memory               ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: R_memory               ! training memory <Rackhun Son>
    dsl4jsb_Real3D_onChunk :: Koppen_memory          ! training memory <Rackhun Son>

    dsl4jsb_Real2D_onChunk :: a_temperature          ! daily average of temperature <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: a_rel_humidity         ! daily average of relative humidity <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: fba                    ! daily fraction of burned area <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: fba_smooth             ! daily smoothed fraction of burned area <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: population             ! daily population density <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: gdp                    ! daily GDP density <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: hdi                    ! daily HDI density <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_primf         ! daily land_use_primf <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_c4per         ! daily land_use_c4per <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_urban         ! daily land_use_urban <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_secdf         ! daily land_use_secdf <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_secma         ! daily land_use_secma <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_primn         ! daily land_use_primn <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_secdn         ! daily land_use_secdn <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_c4ann         ! daily land_use_c4ann <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_pastr         ! daily land_use_pastr <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_range         ! daily land_use_range <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_c3ann         ! daily land_use_c3ann <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_c3nfx         ! daily land_use_c3nfx <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_c3per         ! daily land_use_c3per <Rackhun Son>
    dsl4jsb_Real2D_onChunk :: land_use_secmb         ! daily land_use_secmb <Rackhun Son> 

    dsl4jsb_Real2D_onChunk :: koppen
    dsl4jsb_Real2D_onChunk :: latitude
    dsl4jsb_Real2D_onChunk :: burned_fract
    dsl4jsb_Real2D_onChunk :: burned_fract_infer
    dsl4jsb_Real2D_onChunk :: burned_fract_test
    dsl4jsb_Real2D_onChunk :: damaged_fract

    dsl4jsb_Real2D_onChunk :: co2flux_fire_all_2_atm
    dsl4jsb_Real2D_onChunk :: cflux_dist_green_2_soil
    dsl4jsb_Real2D_onChunk :: cflux_dist_woods_2_soil

    dsl4jsb_Real2D_onChunk :: cover_fract_pot

    !REAL(wp), DIMENSION(options%nc, 9) :: &
    REAL(wp), DIMENSION(options%nc, 13) :: &
      & modified_pfts_fract                ! modified pft fractions for Fire <Rackhun Son>

    ! Locally allocated vectors
    REAL(wp), DIMENSION(options%nc) :: &
      & q_rel_air
         
    ! Declare pointers to variables in memory <Rackhun Son>
    REAL(wp), POINTER :: lightning_climatology(:)
    REAL(wp), POINTER :: slope(:)
    REAL(wp), POINTER :: roughness(:)
    REAL(wp), POINTER :: road_density(:)
    REAL(wp), POINTER :: gfed_region(:)

    REAL(wp), POINTER :: lat(:)

    INTEGER :: iteration, is_training, is_saving, is_infer, memory_idx, n, i_temporal, i_spatial, max_pft_idx
    INTEGER :: save_batch, n_zero, n_nozero
    REAL(wp) :: weight, u, r, t, k, y
    REAL(wp) :: threshold(14)
    CHARACTER(len=10) :: max_pft

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
    iyr     = options%current_datetime%date%year        ! <Rackhun Son> 
    imo     = INT(options%current_datetime%date%month)  ! <Rackhun Son> 
    idy     = options%current_datetime%date%day         ! <Rackhun Son>


    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(DISTURB_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id) ! <Rackhun Son>

    lat => grid%lat(ics:ice, iblk) ! <Rackhun Son>

    ! Get process configurations
    dsl4jsb_Get_config(DISTURB_)

    ! Get process memories
    dsl4jsb_Get_memory(DISTURB_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(FUEL_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(PHENO_)       ! <Rackhun Son>
    dsl4jsb_Get_memory(HYDRO_)       ! <Rackhun Son>
    dsl4jsb_Get_memory(ANTHRO2L_)    ! <Rackhun Son>

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
    dsl4jsb_Get_var2D_onChunk(HYDRO_,     fract_snow)               ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(HYDRO_,     elevation)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(HYDRO_,     wtr_soil_sl)              ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(HYDRO_,     soil_depth_sl)            ! in <Rackhun Son>

    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  a_temperature)            ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  a_rel_humidity)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  fba)                      ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  fba_smooth)               ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  population)               ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  gdp)                      ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  hdi)                      ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_primf)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_c4per)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_urban)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_secdf)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_secma)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_primn)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_secdn)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_c4ann)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_pastr)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_range)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_c3ann)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_c3nfx)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_c3per)           ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(ANTHRO2L_,  land_use_secmb)           ! in <Rackhun Son>
    
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   t_air_ave)                ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_ave)                    ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   q_rel_ave)                ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   ws10_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   rain_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl1_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl2_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl3_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl4_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   swl5_ave)                 ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   lai_ave)                  ! in <Rackhun Son>

    dsl4jsb_Get_var3D_onChunk(DISTURB_,   t_air_annual)             ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   rain_annual)              ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   t_air_clim)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   q_clim)                   ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   q_rel_clim)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl1_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl2_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl3_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl4_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   swl5_clim)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   lai_clim)                 ! in <Rackhun Son>

    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V1_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V2_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V3_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V4_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V5_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V6_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V7_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V8_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V9_memory)                ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V10_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V11_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V12_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V13_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V14_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V15_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V16_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V17_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V18_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V19_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V20_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V21_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V22_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V23_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V24_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V25_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V26_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V27_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V28_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V29_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V30_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V31_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V32_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V33_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V34_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V35_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V36_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V37_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V38_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V39_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V40_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V41_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V42_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V43_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V44_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V45_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V46_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V47_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V48_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V49_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V50_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V51_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V52_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   V53_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   RHC_memory)               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   Y_memory  )               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   R_memory  )               ! in <Rackhun Son>
    dsl4jsb_Get_var3D_onChunk(DISTURB_,   Koppen_memory  )          ! in <Rackhun Son>
    
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   koppen)                   ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   latitude  )               ! in <Rackhun Son>
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   burned_fract)             ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   burned_fract_infer)       ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   burned_fract_test)        ! out
    dsl4jsb_Get_var2D_onChunk(DISTURB_,   damaged_fract)            ! out

    dsl4jsb_Get_var2D_onChunk(CARBON_,    co2flux_fire_all_2_atm)   ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_green_2_soil)  ! out
    dsl4jsb_Get_var2D_onChunk(CARBON_,    cflux_dist_woods_2_soil)  ! out
    
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,    cover_fract_pot)          ! in
    ! -------------------------------------------------------------------------------------------------- !
    !print*, shape(R_memory)


    !q_rel_air(:) = get_relative_humidity_air(nc, q_air(:), t_unfilt(:), press_srf(:))
    !q_rel_air(:) = get_relative_humidity_air_v2(nc, q_air(:), t_air(:), press_srf(:))
    q_rel_air(:) = a_rel_humidity(:)

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
    !DO i = 1, size(soil_depth_sl, dim=1)
    !  IF (soil_depth_sl(i,1) > 0._wp) THEN
    !    swl1(i) = wtr_soil_sl(i,1) / soil_depth_sl(i,1)
    !  ELSE
    !    swl1(i) = 0._wp
    !  ENDIF

    ! calculate daily statistics <Rackhun Son>
    ! get day of year for cliamtology calculation
    day_of_year = INT(get_year_day(options%current_datetime))

    ! skip 2/29 (2/29 is 60th at proleptic gregorian)
    IF (mod(iyr,4) /= 0) THEN
      IF (day_of_year >= 60) day_of_year = day_of_year+1
    ENDIF

    !where (q_rel_air > 100) q_rel_air = 100.

    t_air_ave(:) = t_air_ave(:) + a_temperature(:)*(options%dtime/60._wp/60._wp/24._wp)
    !t_air_ave(:) = t_air_ave(:) + t_air(:)        *(options%dtime/60._wp/60._wp/24._wp)
    q_ave(:)     = q_ave(:)     + q_air(:)        *(options%dtime/60._wp/60._wp/24._wp)
    q_rel_ave(:) = q_rel_ave(:) + q_rel_air(:)    *(options%dtime/60._wp/60._wp/24._wp)
    ws10_ave(:)  = ws10_ave(:)  + wind_10m(:)     *(options%dtime/60._wp/60._wp/24._wp)
    rain_ave(:)  = rain_ave(:)  + rain(:)         *(options%dtime/60._wp/60._wp/24._wp)
    swl1_ave(:)  = swl1_ave(:)  + wtr_soil_sl(:,1)*(options%dtime/60._wp/60._wp/24._wp)
    swl2_ave(:)  = swl2_ave(:)  + wtr_soil_sl(:,2)*(options%dtime/60._wp/60._wp/24._wp)
    swl3_ave(:)  = swl3_ave(:)  + wtr_soil_sl(:,3)*(options%dtime/60._wp/60._wp/24._wp)
    swl4_ave(:)  = swl4_ave(:)  + wtr_soil_sl(:,4)*(options%dtime/60._wp/60._wp/24._wp)
    swl5_ave(:)  = swl5_ave(:)  + wtr_soil_sl(:,5)*(options%dtime/60._wp/60._wp/24._wp)
    lai_ave(:)   = lai_ave(:)   + lai(:)          *(options%dtime/60._wp/60._wp/24._wp)

    IF (is_newday(options%current_datetime, options%dtime)) THEN ! only once per day
      threshold(1)=0.;   threshold(2)=-0.5;  threshold(3)=-0.5;  threshold(4)=0.;  threshold(5)=-0.5 
      threshold(6)=0.;   threshold(7)=0.;  threshold(8)=0.75;  threshold(9)=0.75;  threshold(10)=0. 
      threshold(11)=0.;  threshold(12)=0.; threshold(13)=0.; threshold(14)=0.

      ! calculate auxiliary variables for DL fire model   <Rackhun Son>
      IF (dsl4jsb_Config(DISTURB_)%fire_algorithm == 5) THEN
        diff_year = FLOAT(iyr-time_config%tc_exp_startdate%date%year) +1._wp

        iteration  = 0
        IF ((iyr .ge. 2003) .AND. (iyr .le. 2012)) THEN
        !IF (iyr .eq. 2003) THEN
          IF (iteration .eq. 0) THEN
            CALL get_first_iblk(iblk)
          ENDIF

          is_training = 1
          !iteration = (iyr-2003)*5946*365 + (day_of_year-1)*5946 + iblk*600
          iteration = (iyr-2003)*365 + day_of_year
        ELSE
          is_training = 0
        ENDIF

        IF ((iyr .ge. 2001) .AND. (iyr .le. 2012)) THEN
        !IF ((iyr .ge. 2001) .AND. (iyr .le. 2005)) THEN
          is_saving = 1
          memory_idx = day_of_year + MOD(iyr,2)*365

          IF (memory_idx .gt. 730) THEN
            is_saving = 0
          ENDIF

          !max_pft_idx = MAXLOC(pfts_fract(i_spatial,:), dim=1)
          max_pft_idx = MAXLOC(SUM(pfts_fract(:,:), dim=1), dim=1)
          IF (max_pft_idx .eq. 1) THEN
            max_pft = 'pft01'
          ELSEIF (max_pft_idx .eq. 2) THEN
            max_pft = 'pft02'
          ELSEIF (max_pft_idx .eq. 3) THEN
            max_pft = 'pft03'
          ELSEIF (max_pft_idx .eq. 4) THEN
            max_pft = 'pft04'
          ELSEIF (max_pft_idx .eq. 5) THEN
            max_pft = 'pft05'
          ELSEIF (max_pft_idx .eq. 6) THEN
            max_pft = 'pft06'
          ELSEIF (max_pft_idx .eq. 7) THEN
            max_pft = 'pft07'
          ELSEIF (max_pft_idx .eq. 8) THEN
            max_pft = 'pft08'
          ELSEIF (max_pft_idx .eq. 9) THEN
            max_pft = 'pft09'
          ELSEIF (max_pft_idx .eq. 10) THEN
            max_pft = 'pft10'
          ELSEIF (max_pft_idx .eq. 11) THEN
            max_pft = 'pft11'
          ENDIF

          IF  (iyr .ge. 2002) THEN
            ! random sampling from memory
            r = SUM(R_memory)
            y = SUM(Y_memory)
            !k = SUM(Koppen_memory)
            IF (tile%name .eq. max_pft) THEN
              !IF ((r .gt. 0.) .and. (y .gt. 0.)) THEN
            !IF ((tile%name .eq. max_pft) .and. (r .gt. 0.) .and. (k .gt. 0.)) THEN
                
              n=1
              n_zero=0
              n_nozero=0
              DO WHILE ( n<=8 )
                CALL RANDOM_NUMBER(u)
                i_temporal = 1 + INT(FLOOR((730)*u))
                CALL RANDOM_NUMBER(u)
                i_spatial  = 1 + INT(FLOOR((128)*u))

                IF (R_memory(i_spatial,i_temporal) .ge. 1) THEN
                  IF (Y_memory(i_spatial,i_temporal) .gt. 0) THEN
                    n_nozero = n_nozero+1
                    !IF (n_nozero .le. 8) THEN
                    save_batch = 1
                    !ELSE
                    !  save_batch = 0
                    !ENDIF
                  ELSE
                    CALL RANDOM_NUMBER(u)
                    IF (u .gt. 0.5) THEN
                      save_batch = 1
                    ELSE
                      save_batch = 0
                    ENDIF
                  ENDIF
                  !ELSE
                  !  n_zero = n_zero+1
                  !  IF (n_zero .le. 4) THEN
                  !    save_batch = 1
                  !  ELSE
                  !    save_batch = 0
                  !  ENDIF
                  !ENDIF
                  !save_batch = 1

                  IF (save_batch .eq. 1) THEN
                    CALL update_pid_batch(                                                                                    &
                      & V1_memory (i_spatial,i_temporal),  V2_memory(i_spatial,i_temporal),  V3_memory(i_spatial,i_temporal), &
                      & V4_memory (i_spatial,i_temporal),  V5_memory(i_spatial,i_temporal),  V6_memory(i_spatial,i_temporal), &
                      & V7_memory (i_spatial,i_temporal),  V8_memory(i_spatial,i_temporal),  V9_memory(i_spatial,i_temporal), &
                      & V10_memory(i_spatial,i_temporal), V11_memory(i_spatial,i_temporal), V12_memory(i_spatial,i_temporal), &
                      & V13_memory(i_spatial,i_temporal), V14_memory(i_spatial,i_temporal), V15_memory(i_spatial,i_temporal), &
                      & V16_memory(i_spatial,i_temporal), V17_memory(i_spatial,i_temporal), V18_memory(i_spatial,i_temporal), &
                      & V19_memory(i_spatial,i_temporal), V20_memory(i_spatial,i_temporal), V21_memory(i_spatial,i_temporal), &
                      & V22_memory(i_spatial,i_temporal), V23_memory(i_spatial,i_temporal), V24_memory(i_spatial,i_temporal), &
                      & V25_memory(i_spatial,i_temporal), V26_memory(i_spatial,i_temporal), V27_memory(i_spatial,i_temporal), &
                      & V28_memory(i_spatial,i_temporal), V29_memory(i_spatial,i_temporal), V30_memory(i_spatial,i_temporal), &
                      & V31_memory(i_spatial,i_temporal), V32_memory(i_spatial,i_temporal), V33_memory(i_spatial,i_temporal), &
                      & V34_memory(i_spatial,i_temporal), V35_memory(i_spatial,i_temporal), V36_memory(i_spatial,i_temporal), &
                      & V37_memory(i_spatial,i_temporal), V38_memory(i_spatial,i_temporal), V39_memory(i_spatial,i_temporal), &
                      & V40_memory(i_spatial,i_temporal), V41_memory(i_spatial,i_temporal), V42_memory(i_spatial,i_temporal), &
                      & V43_memory(i_spatial,i_temporal), V44_memory(i_spatial,i_temporal), V45_memory(i_spatial,i_temporal), &
                      & V46_memory(i_spatial,i_temporal), V47_memory(i_spatial,i_temporal), V48_memory(i_spatial,i_temporal), &
                      & V49_memory(i_spatial,i_temporal), V50_memory(i_spatial,i_temporal), V51_memory(i_spatial,i_temporal), &
                      & V52_memory(i_spatial,i_temporal), V53_memory(i_spatial,i_temporal),                                   &
                      & RHC_memory(i_spatial,i_temporal), Y_memory  (i_spatial,i_temporal),   R_memory(i_spatial,i_temporal), &
                      & Koppen_memory(i_spatial,i_temporal), n, iblk)
                    n = n+1
                  ENDIF
                  !ENDDO

                ELSE
                  !DO n = 1,8
                  CALL update_pid_batch(                                                                                    &
                    & 0._wp,  0._wp,  0._wp, & ! 1
                    & 0._wp,  0._wp,  0._wp, & ! 2
                    & 0._wp,  0._wp,  0._wp, & ! 3
                    & 0._wp,  0._wp,  0._wp, & ! 4
                    & 0._wp,  0._wp,  0._wp, & ! 5
                    & 0._wp,  0._wp,  0._wp, & ! 6
                    & 0._wp,  0._wp,  0._wp, & ! 7
                    & 0._wp,  0._wp,  0._wp, & ! 8
                    & 0._wp,  0._wp,  0._wp, & ! 9
                    & 0._wp,  0._wp,  0._wp, & ! 10
                    & 0._wp,  0._wp,  0._wp, & ! 11
                    & 0._wp,  0._wp,  0._wp, & ! 12
                    & 0._wp,  0._wp,  0._wp, & ! 13
                    & 0._wp,  0._wp,  0._wp, & ! 14
                    & 0._wp,  0._wp,  0._wp, & ! 15
                    & 0._wp,  0._wp,  0._wp, & ! 16
                    & 0._wp,  0._wp,  0._wp, & ! 17
                    & 0._wp,  0._wp,         & ! 18

                    & 0._wp,  0._wp,  0._wp, & ! 19
                    & 0._wp, n, iblk)
                  n = n+1
                  !ENDDO
                ENDIF
                !save_batch = 1
                !t = threshold(INT(R_memory(i_spatial,i_temporal)))
                !IF (t .ge. 0) THEN
                !  IF (Y_memory(i_spatial,i_temporal) .gt. 0.) THEN
                !    save_batch = 1
                !  ELSE
                !    CALL RANDOM_NUMBER(u)
                !    IF (u .gt. t) THEN
                !      save_batch = 1
                !    ELSE
                !      save_batch = 0
                !    ENDIF
                !  ENDIF
                !ELSE
                !  IF (Y_memory(i_spatial,i_temporal) .eq. 0.) THEN
                !    save_batch = 1
                !  ELSE
                !    CALL RANDOM_NUMBER(u)
                !    IF (u .gt. -t) THEN
                !      save_batch = 1
                !    ELSE
                !      save_batch = 0
                !    ENDIF
                !  ENDIF
                !ENDIF

                !IF (save_batch .eq. 1) THEN
                !  CALL update_pid_batch(                                                                                    &
                !    & V1_memory (i_spatial,i_temporal),  V2_memory(i_spatial,i_temporal),  V3_memory(i_spatial,i_temporal), &
                !    & V4_memory (i_spatial,i_temporal),  V5_memory(i_spatial,i_temporal),  V6_memory(i_spatial,i_temporal), &
                !    & V7_memory (i_spatial,i_temporal),  V8_memory(i_spatial,i_temporal),  V9_memory(i_spatial,i_temporal), &
                !    & V10_memory(i_spatial,i_temporal), V11_memory(i_spatial,i_temporal), V12_memory(i_spatial,i_temporal), &
                !    & V13_memory(i_spatial,i_temporal), V14_memory(i_spatial,i_temporal), V15_memory(i_spatial,i_temporal), &
                !    & V16_memory(i_spatial,i_temporal), V17_memory(i_spatial,i_temporal), V18_memory(i_spatial,i_temporal), &
                !    & V19_memory(i_spatial,i_temporal), V20_memory(i_spatial,i_temporal), V21_memory(i_spatial,i_temporal), &
                !    & V22_memory(i_spatial,i_temporal), V23_memory(i_spatial,i_temporal), V24_memory(i_spatial,i_temporal), &
                !    & V25_memory(i_spatial,i_temporal), V26_memory(i_spatial,i_temporal), V27_memory(i_spatial,i_temporal), &
                !    & V28_memory(i_spatial,i_temporal), V29_memory(i_spatial,i_temporal), V30_memory(i_spatial,i_temporal), &
                !    & V31_memory(i_spatial,i_temporal), V32_memory(i_spatial,i_temporal), V33_memory(i_spatial,i_temporal), &
                !    & V34_memory(i_spatial,i_temporal), V35_memory(i_spatial,i_temporal), V36_memory(i_spatial,i_temporal), &
                !    & V37_memory(i_spatial,i_temporal), V38_memory(i_spatial,i_temporal), V39_memory(i_spatial,i_temporal), &
                !    & RHC_memory(i_spatial,i_temporal), Y_memory  (i_spatial,i_temporal),   R_memory(i_spatial,i_temporal), &
                !    & Koppen_memory(i_spatial,i_temporal), n, iblk)
                !  n = n+1
                !ENDIF
              ENDDO
            ENDIF
          ENDIF
          
        ELSE
          is_saving = 0
        ENDIF

        IF (iyr .ge. 2013) THEN
        !IF (iyr .ge. 2005) THEN
          is_infer = 1
        ELSE
          is_infer = 0
        ENDIF

        !IF ((iyr .ge. 2002) .AND. (iyr .le. 2012)) THEN
        !IF ((iyr .ge. 2003) .AND. (iyr .le. 2005)) THEN

          ! random sampling from memory
          !r = SUM(R_memory)
          !IF ((tile%name .eq. max_pft) .and. (r .gt. 0.)) THEN
          !  n=1
          !  DO WHILE ( n<=8 )
          !    CALL RANDOM_NUMBER(u)
          !    i_temporal = 1 + INT(FLOOR((730)*u))
          !    CALL RANDOM_NUMBER(u)
          !    i_spatial  = 1 + INT(FLOOR((128)*u))

              !t = threshold(INT(R_memory(i_spatial,i_temporal)))
              !IF (t .ge. 0) THEN
              !  IF (Y_memory(i_spatial,i_temporal) .gt. 0.) THEN
              !    save_batch = 1
              !  ELSE
              !    CALL RANDOM_NUMBER(u)
              !    IF (u .gt. t) THEN
              !      save_batch = 1
              !    ELSE
              !      save_batch = 0
              !    ENDIF
              !  ENDIF
              !ELSE
              !  IF (Y_memory(i_spatial,i_temporal) .eq. 0.) THEN
              !    save_batch = 1
              !  ELSE
              !    CALL RANDOM_NUMBER(u)
              !    IF (u .gt. -t) THEN
              !      save_batch = 1
              !    ELSE
              !      save_batch = 0
              !    ENDIF
              !  ENDIF
              !ENDIF
              
              !IF ((R_memory(i_spatial,i_temporal) .ne. 0.) .and. (save_batch .eq. 1)) THEN
              !IF ((R_memory(i_spatial,i_temporal) .ne. 0.) .and. (Koppen_memory(i_spatial,i_temporal) .ne. 0.))  THEN
              !  CALL update_pid_batch(                                                                                    &
              !    & V1_memory (i_spatial,i_temporal),  V2_memory(i_spatial,i_temporal),  V3_memory(i_spatial,i_temporal), &
              !    & V4_memory (i_spatial,i_temporal),  V5_memory(i_spatial,i_temporal),  V6_memory(i_spatial,i_temporal), &
              !    & V7_memory (i_spatial,i_temporal),  V8_memory(i_spatial,i_temporal),  V9_memory(i_spatial,i_temporal), &
              !    & V10_memory(i_spatial,i_temporal), V11_memory(i_spatial,i_temporal), V12_memory(i_spatial,i_temporal), &
              !    & V13_memory(i_spatial,i_temporal), V14_memory(i_spatial,i_temporal), V15_memory(i_spatial,i_temporal), &
              !    & V16_memory(i_spatial,i_temporal), V17_memory(i_spatial,i_temporal), V18_memory(i_spatial,i_temporal), &
              !    & V19_memory(i_spatial,i_temporal), V20_memory(i_spatial,i_temporal), V21_memory(i_spatial,i_temporal), &
              !    & V22_memory(i_spatial,i_temporal), V23_memory(i_spatial,i_temporal), V24_memory(i_spatial,i_temporal), &
              !    & V25_memory(i_spatial,i_temporal), V26_memory(i_spatial,i_temporal), V27_memory(i_spatial,i_temporal), &
              !    & V28_memory(i_spatial,i_temporal), V29_memory(i_spatial,i_temporal), V30_memory(i_spatial,i_temporal), &
              !    & V31_memory(i_spatial,i_temporal), V32_memory(i_spatial,i_temporal), V33_memory(i_spatial,i_temporal), &
              !    & V34_memory(i_spatial,i_temporal), V35_memory(i_spatial,i_temporal), V36_memory(i_spatial,i_temporal), &
              !    & V37_memory(i_spatial,i_temporal), V38_memory(i_spatial,i_temporal), V39_memory(i_spatial,i_temporal), &
              !    & RHC_memory(i_spatial,i_temporal), Y_memory  (i_spatial,i_temporal),   R_memory(i_spatial,i_temporal), &
              !    & Koppen_memory(i_spatial,i_temporal), n, iblk)
              !  n = n+1
              !ENDIF
          !  ENDDO
          !ENDIF
        !ENDIF

        IF (is_time_experiment_start(options%current_datetime)) THEN   ! 1st iteration
          t_air_ave(:) = t_air_ave(:) /(options%dtime/60._wp/60._wp/24._wp)
          q_ave(:)     = q_ave(:)     /(options%dtime/60._wp/60._wp/24._wp)
          q_rel_ave(:) = q_rel_ave(:) /(options%dtime/60._wp/60._wp/24._wp)
          ws10_ave(:)  = ws10_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          rain_ave(:)  = rain_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          swl1_ave(:)  = swl1_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          swl2_ave(:)  = swl2_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          swl3_ave(:)  = swl3_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          swl4_ave(:)  = swl4_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          swl5_ave(:)  = swl5_ave(:)  /(options%dtime/60._wp/60._wp/24._wp)
          lai_ave(:)   = lai_ave(:)   /(options%dtime/60._wp/60._wp/24._wp)

          t_air_clim(:,day_of_year) = t_air_ave(:) !/(options%dtime/60._wp/60._wp/24._wp)
          q_clim(:,day_of_year)     = q_ave(:)     !/(options%dtime/60._wp/60._wp/24._wp)
          q_rel_clim(:,day_of_year) = q_rel_ave(:) !/(options%dtime/60._wp/60._wp/24._wp)
          swl1_clim(:,day_of_year)  = swl1_ave(:)  !/(options%dtime/60._wp/60._wp/24._wp)
          swl2_clim(:,day_of_year)  = swl2_ave(:)  !/(options%dtime/60._wp/60._wp/24._wp)
          swl3_clim(:,day_of_year)  = swl3_ave(:)  !/(options%dtime/60._wp/60._wp/24._wp)
          swl4_clim(:,day_of_year)  = swl4_ave(:)  !/(options%dtime/60._wp/60._wp/24._wp)
          swl5_clim(:,day_of_year)  = swl5_ave(:)  !/(options%dtime/60._wp/60._wp/24._wp)
          lai_clim(:,day_of_year)   = lai_ave(:)   !/(options%dtime/60._wp/60._wp/24._wp)


          ! initialize trainable parameters
          CALL update_trainable_parameters(                             &
            & dsl4jsb_Config(DISTURB_)%layer1_weight,                   &
            & dsl4jsb_Config(DISTURB_)%layer1_bias,                     &
            & dsl4jsb_Config(DISTURB_)%layer2_weight,                   &
            & dsl4jsb_Config(DISTURB_)%layer2_bias,                     &
            & dsl4jsb_Config(DISTURB_)%layer3_weight,                   &
            & dsl4jsb_Config(DISTURB_)%layer3_bias,                     &
            & dsl4jsb_Config(DISTURB_)%layer4_weight,                   &
            & dsl4jsb_Config(DISTURB_)%layer4_bias                      )

        !ELSEIF (iyr == time_config%tc_exp_startdate%date%year) THEN    ! 1st year
        ELSE

          IF (iyr == time_config%tc_exp_startdate%date%year) THEN    ! 1st year
            t_air_clim(:,day_of_year) = t_air_ave(:)
            q_clim(:,day_of_year)     = q_ave(:)
            q_rel_clim(:,day_of_year) = q_rel_ave(:)
            swl1_clim(:,day_of_year)  = swl1_ave(:)
            swl2_clim(:,day_of_year)  = swl2_ave(:)
            swl3_clim(:,day_of_year)  = swl3_ave(:)
            swl4_clim(:,day_of_year)  = swl4_ave(:)
            swl5_clim(:,day_of_year)  = swl5_ave(:)
            lai_clim(:,day_of_year)   = lai_ave(:)

          ELSEIF (day_of_year == 60) Then                                ! 2/29 <= 2/28 (always)
            t_air_clim(:,day_of_year) = t_air_clim(:,59)
            q_clim(:,day_of_year)     = q_clim(:,59)
            q_rel_clim(:,day_of_year) = q_rel_clim(:,59)
            swl1_clim(:,day_of_year)  = swl1_clim(:,59)
            swl2_clim(:,day_of_year)  = swl2_clim(:,59)
            swl3_clim(:,day_of_year)  = swl3_clim(:,59)
            swl4_clim(:,day_of_year)  = swl4_clim(:,59)
            lai_clim(:,day_of_year)   = lai_clim(:,59)

          ELSE
            t_air_clim(:,day_of_year) = t_air_clim(:,day_of_year) + t_air_ave(:)
            q_clim(:,day_of_year)     = q_clim(:,day_of_year)     + q_ave(:)
            q_rel_clim(:,day_of_year) = q_rel_clim(:,day_of_year) + q_rel_ave(:)
            swl1_clim(:,day_of_year)  = swl1_clim(:,day_of_year)  + swl1_ave(:)
            swl2_clim(:,day_of_year)  = swl2_clim(:,day_of_year)  + swl2_ave(:)
            swl3_clim(:,day_of_year)  = swl3_clim(:,day_of_year)  + swl3_ave(:)
            swl4_clim(:,day_of_year)  = swl4_clim(:,day_of_year)  + swl4_ave(:)
            swl5_clim(:,day_of_year)  = swl5_clim(:,day_of_year)  + swl5_ave(:)
            lai_clim(:,day_of_year)   = lai_clim(:,day_of_year)   + lai_ave(:)
          ENDIF

          IF ((iyr .eq. 2013) .AND. (imo .eq. 1) .AND. (idy .eq. 1)) THEN
          !IF ((iyr .eq. 2005) .AND. (imo .eq. 1) .AND. (idy .eq. 1)) THEN
            CALL load_trained_parameters(                               &
              & dsl4jsb_Config(DISTURB_)%layer1_weight_update,          &
              & dsl4jsb_Config(DISTURB_)%layer1_bias_update,            &
              & dsl4jsb_Config(DISTURB_)%layer2_weight_update,          &
              & dsl4jsb_Config(DISTURB_)%layer2_bias_update,            &
              & dsl4jsb_Config(DISTURB_)%layer3_weight_update,          &
              & dsl4jsb_Config(DISTURB_)%layer3_bias_update,            &
              & dsl4jsb_Config(DISTURB_)%layer4_weight_update,          &
              & dsl4jsb_Config(DISTURB_)%layer4_bias_update             )
            print*, dsl4jsb_Config(DISTURB_)%layer1_weight_update(10,128*2+20)
          ENDIF

        ENDIF

        ! modify pft fractions for Fire
        modified_pfts_fract(:,1)  = fract_snow(:)
        modified_pfts_fract(:,2)  = pfts_fract(:,1) *(1-fract_snow(:))*fract_fpc_max     ! Tropical broadleaf evergreen
        modified_pfts_fract(:,3)  = pfts_fract(:,2) *(1-fract_snow(:))*fract_fpc_max     ! Tropical broadleaf deciduous
        modified_pfts_fract(:,4)  = pfts_fract(:,3) *(1-fract_snow(:))*fract_fpc_max     ! Extra-tropical broadleaf evergreen
        modified_pfts_fract(:,5)  = pfts_fract(:,4) *(1-fract_snow(:))*fract_fpc_max     ! Extra_tropical broadleaf deciduous
        modified_pfts_fract(:,6)  = pfts_fract(:,5) *(1-fract_snow(:))*fract_fpc_max     ! Raingreen Shrubs
        modified_pfts_fract(:,7)  = pfts_fract(:,6) *(1-fract_snow(:))*fract_fpc_max     ! Deciduous Shrubs
        !modified_pfts_fract(:,8)  = SUM(pfts_fract(:,7:8),DIM=2)*(1-fract_snow(:))*fract_fpc_max   ! Grass
        modified_pfts_fract(:,8)  = pfts_fract(:,7) *(1-fract_snow(:))*fract_fpc_max     !
        modified_pfts_fract(:,9)  = pfts_fract(:,8) *(1-fract_snow(:))*fract_fpc_max     !
        modified_pfts_fract(:,10) = pfts_fract(:,9) *(1-fract_snow(:))*fract_fpc_max     !
        modified_pfts_fract(:,11) = pfts_fract(:,10)*(1-fract_snow(:))*fract_fpc_max     !
        modified_pfts_fract(:,12) = pfts_fract(:,11)*(1-fract_snow(:))*fract_fpc_max     !
        modified_pfts_fract(:,13) = (1-fract_snow(:))*(1-fract_fpc_max)                  ! Others

        lightning_climatology => dsl4jsb_memory(DISTURB_)%lightning_climatology(ics:ice, iblk, day_of_year)    ! <Rackhun Son>

        slope                     => dsl4jsb_memory(DISTURB_)%slope(ics:ice, iblk)                     ! <Rackhun Son>
        roughness                 => dsl4jsb_memory(DISTURB_)%roughness(ics:ice, iblk)                 ! <Rackhun Son>
        road_density              => dsl4jsb_memory(DISTURB_)%road_density(ics:ice, iblk)              ! <Rackhun Son>
        gfed_region               => dsl4jsb_memory(DISTURB_)%gfed_region(ics:ice, iblk)               ! <Rackhun Son>

        ! annual record & Koppen_Geiger_classification
        t_air_annual(:,day_of_year) = t_air_ave(:)
        rain_annual (:,day_of_year) = rain_ave(:)

        CALL Koppen_Geiger_classification_lv2(nc, t_air_annual, rain_annual, lat, koppen)
        latitude(:) = lat(:)

      ENDIF

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

            CALL burned_fract_dl(                                         &
              & t_air_ave(:),                                             &
              & t_air_clim(:,day_of_year)/diff_year,                      &  
              & q_rel_ave(:),                                             &
              & q_rel_clim(:,day_of_year)/diff_year,                      &
              & q_ave(:),                                                 &
              & q_clim(:,day_of_year)/diff_year,                          &
              & ws10_ave(:),                                              &
              & rain_ave(:),                                              &
              & swl1_ave(:),                                              &
              & swl1_clim(:,day_of_year)/diff_year,                       &
              & swl2_ave(:),                                              &
              & swl2_clim(:,day_of_year)/diff_year,                       &
              & swl3_ave(:),                                              &
              & swl3_clim(:,day_of_year)/diff_year,                       &
              & swl4_ave(:),                                              &
              & swl4_clim(:,day_of_year)/diff_year,                       &
              & swl5_ave(:),                                              &
              & swl5_clim(:,day_of_year)/diff_year,                       &
              & lai_ave(:),                                               &
              & lai_clim(:,day_of_year)/diff_year,                        &
              & modified_pfts_fract(:,:),                                 &
              !& pfts_fract(:,:),                                          &
              & fuel(:),                                                  &
              & lightning_climatology(:),                                 &
              & elevation(:),                                             &
              !& slope(:),                                                 &
              !& roughness(:),                                             &

              & population(:),                                            &
              & gdp(:),                                                   &
              & hdi(:),                                                   &

              & land_use_c3ann(:),                                        & ! in
              & land_use_c3nfx(:),                                        & ! in
              & land_use_c3per(:),                                        & ! in
              & land_use_c4ann(:),                                        & ! in
              & land_use_c4per(:),                                        & ! in
              & land_use_pastr(:),                                        & ! in
              & land_use_primf(:),                                        & ! in
              & land_use_primn(:),                                        & ! in
              & land_use_range(:),                                        & ! in
              & land_use_secdf(:),                                        & ! in
              & land_use_secdn(:),                                        & ! in
              & land_use_secma(:),                                        & ! in
              & land_use_secmb(:),                                        & ! in
              & land_use_urban(:),                                        & ! in

              !& t_air_annual(:,:),                                        &
              !& rain_annual(:,:),                                         &
              & koppen(:),                                                &

              !& fba,                                                      &
              & fba_smooth,                                               &
              & nc,                                                       &
              & iblk,                                                     &
              & tile%name,                                                &
              & max_pft,                                                  &
              
              & dsl4jsb_Config(DISTURB_)%d_input,                         &
              & dsl4jsb_Config(DISTURB_)%d_hidden,                        &
              & dsl4jsb_Config(DISTURB_)%d_output,                        &
              & dsl4jsb_Config(DISTURB_)%d_region,                        &
              & dsl4jsb_Config(DISTURB_)%d_class,                         &

              & dsl4jsb_Config(DISTURB_)%layer1_weight,                   &
              & dsl4jsb_Config(DISTURB_)%layer1_bias,                     &
              & dsl4jsb_Config(DISTURB_)%layer2_weight,                   &
              & dsl4jsb_Config(DISTURB_)%layer2_bias,                     &
              & dsl4jsb_Config(DISTURB_)%layer3_weight,                   &
              & dsl4jsb_Config(DISTURB_)%layer3_bias,                     &
              & dsl4jsb_Config(DISTURB_)%layer4_weight,                   &
              & dsl4jsb_Config(DISTURB_)%layer4_bias,                     &

              & dsl4jsb_Config(DISTURB_)%layer1_weight_update,            &
              & dsl4jsb_Config(DISTURB_)%layer1_bias_update,              &
              & dsl4jsb_Config(DISTURB_)%layer2_weight_update,            &
              & dsl4jsb_Config(DISTURB_)%layer2_bias_update,              &
              & dsl4jsb_Config(DISTURB_)%layer3_weight_update,            &
              & dsl4jsb_Config(DISTURB_)%layer3_bias_update,              &
              & dsl4jsb_Config(DISTURB_)%layer4_weight_update,            &
              & dsl4jsb_Config(DISTURB_)%layer4_bias_update,              &

              & dsl4jsb_Config(DISTURB_)%mW1,                             &
              & dsl4jsb_Config(DISTURB_)%mb1,                             &
              & dsl4jsb_Config(DISTURB_)%mW2,                             &
              & dsl4jsb_Config(DISTURB_)%mb2,                             &
              & dsl4jsb_Config(DISTURB_)%mW3,                             &
              & dsl4jsb_Config(DISTURB_)%mb3,                             &
              & dsl4jsb_Config(DISTURB_)%mW4,                             &
              & dsl4jsb_Config(DISTURB_)%mb4,                             &

              & dsl4jsb_Config(DISTURB_)%vW1,                             &
              & dsl4jsb_Config(DISTURB_)%vb1,                             &
              & dsl4jsb_Config(DISTURB_)%vW2,                             &
              & dsl4jsb_Config(DISTURB_)%vb2,                             &
              & dsl4jsb_Config(DISTURB_)%vW3,                             &
              & dsl4jsb_Config(DISTURB_)%vb3,                             &
              & dsl4jsb_Config(DISTURB_)%vW4,                             &
              & dsl4jsb_Config(DISTURB_)%vb4,                             &

              & iteration,                                                &
              & is_training,                                              &
              & is_saving,                                                &
              & is_infer,                                                 &
              & memory_idx,                                               &

              & V1_memory (:,:),  V2_memory(:,:),  V3_memory(:,:),        &
              & V4_memory (:,:),  V5_memory(:,:),  V6_memory(:,:),        &
              & V7_memory (:,:),  V8_memory(:,:),  V9_memory(:,:),        &
              & V10_memory(:,:), V11_memory(:,:), V12_memory(:,:),        &
              & V13_memory(:,:), V14_memory(:,:), V15_memory(:,:),        &
              & V16_memory(:,:), V17_memory(:,:), V18_memory(:,:),        &
              & V19_memory(:,:), V20_memory(:,:), V21_memory(:,:),        &
              & V22_memory(:,:), V23_memory(:,:), V24_memory(:,:),        &
              & V25_memory(:,:), V26_memory(:,:), V27_memory(:,:),        &
              & V28_memory(:,:), V29_memory(:,:), V30_memory(:,:),        &
              & V31_memory(:,:), V32_memory(:,:), V33_memory(:,:),        &
              & V34_memory(:,:), V35_memory(:,:), V36_memory(:,:),        &
              & V37_memory(:,:), V38_memory(:,:), V39_memory(:,:),        &
              & V40_memory(:,:), V41_memory(:,:), V42_memory(:,:),        &
              & V43_memory(:,:), V44_memory(:,:), V45_memory(:,:),        &
              & V46_memory(:,:), V47_memory(:,:), V48_memory(:,:),        &
              & V49_memory(:,:), V50_memory(:,:), V51_memory(:,:),        &
              & V52_memory(:,:), V53_memory(:,:),                         &
              & RHC_memory(:,:),   Y_memory(:,:),   R_memory(:,:),        &
              & Koppen_memory(:,:),                                       &

              & dsl4jsb_Config(DISTURB_)%min_max,                         &
              & gfed_region(:),                                           &

              !& koppen(:),                                                & ! inout
              & burned_fract(:),                                          & ! inout
              & burned_fract_infer(:),                                    & ! inout
              & burned_fract_test(:),                                     & ! inout
              & weight)


              ! save training results
              IF (iyr == 2015) then
              !IF (iyr == 2005) then

                open(1001,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/W1.bin', &
                  access='direct', form='unformatted', &
                  recl=dsl4jsb_Config(DISTURB_)%d_input*dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                open(1002,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/W2.bin', &
                  access='direct', form='unformatted', &
                  recl=dsl4jsb_Config(DISTURB_)%d_hidden*dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                open(1003,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/W3.bin', &
                  access='direct', form='unformatted', &
                  recl=dsl4jsb_Config(DISTURB_)%d_hidden*dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                open(1004,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/W4.bin', &
                  access='direct', form='unformatted', &
                  recl=dsl4jsb_Config(DISTURB_)%d_hidden*dsl4jsb_Config(DISTURB_)%d_output_n_class*4)
                !open(1004,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/W4.bin', &
                !  access='direct', form='unformatted', &
                !  recl=dsl4jsb_Config(DISTURB_)%d_hidden*dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                !open(1005,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/W5.bin', &
                !  access='direct', form='unformatted', &
                !  recl=dsl4jsb_Config(DISTURB_)%d_hidden*dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                !open(1006,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/W6.bin', &
                !  access='direct', form='unformatted', &
                !  recl=dsl4jsb_Config(DISTURB_)%d_hidden*dsl4jsb_Config(DISTURB_)%d_output_n_class*4)

                open(2001,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/b1.bin', &
                  access='direct', form='unformatted', &
                  recl=dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                open(2002,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/b2.bin', &
                  access='direct', form='unformatted', &
                  recl=dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                open(2003,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/b3.bin', &
                  access='direct', form='unformatted', &
                  recl=dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                open(2004,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/b4.bin', &
                  access='direct', form='unformatted', &
                  recl=dsl4jsb_Config(DISTURB_)%d_output_n_class*4)
                !open(2004,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/b4.bin', &
                !  access='direct', form='unformatted', &
                !  recl=dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                !open(2005,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/b5.bin', &
                !  access='direct', form='unformatted', &
                !  recl=dsl4jsb_Config(DISTURB_)%d_hidden_n_class*4)
                !open(2006,file='/work/mj0143/b301058/fire_online_2_koppen_lv2/saved_parameters/b6.bin', &
                !  access='direct', form='unformatted', &
                !  recl=dsl4jsb_Config(DISTURB_)%d_output_n_class*4)

                write(1001,rec=1) layer1_weight_global
                write(1002,rec=1) layer2_weight_global
                write(1003,rec=1) layer3_weight_global
                write(1004,rec=1) layer4_weight_global
                !write(1005,rec=1) layer5_weight_global
                !write(1006,rec=1) layer6_weight_global

                write(2001,rec=1) layer1_bias_global
                write(2002,rec=1) layer2_bias_global
                write(2003,rec=1) layer3_bias_global
                write(2004,rec=1) layer4_bias_global
                !write(2005,rec=1) layer5_bias_global
                !write(2006,rec=1) layer6_bias_global

                print*, layer1_weight_global(10,128*2+20)

              ENDIF


          CASE DEFAULT
            CALL finish('disturbed_fract','Unknown fire algorithm')
        END SELECT
      END IF

      ! reset daily statistics
      t_air_ave(:) = 0.
      q_ave(:) = 0.
      q_rel_ave(:) = 0.
      ws10_ave(:) = 0.
      rain_ave(:) = 0.
      swl1_ave(:) = 0.
      swl2_ave(:) = 0.
      swl3_ave(:) = 0.
      swl4_ave(:) = 0.
      swl5_ave(:) = 0.
      lai_ave(:) = 0.


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
    dsl4jsb_Aggregate_onChunk(DISTURB_, burned_fract_infer,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, burned_fract_test,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, koppen,                  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(DISTURB_, latitude,                weighted_by_fract)
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

    USE mo_carbon_process,        ONLY: relocate_carbon_fire, relocate_carbon_damage
    USE mo_carbon_interface,      ONLY: recalc_per_tile_vars, calculate_current_c_ag_1_and_bg_sums,  &
                                      & calculate_current_c_ta_state_sum, check_carbon_conservation, &
                                      & recalc_carbon_per_tile_vars
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_time,              ONLY: is_newday

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    REAL(wp), DIMENSION(options%nc)           :: old_c_state_sum_ta, dummy_flux

    INTEGER :: iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_natural_disturbances'

    ! Declare process memories
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_config(DISTURB_)
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

    IF ( tile%Is_process_active(CARBON_) ) THEN
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
