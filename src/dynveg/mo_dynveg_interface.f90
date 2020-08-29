!>
!! @brief Contains the interfaces to the dynveg process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2018-12-20)
!!
MODULE mo_dynveg_interface
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
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_io_netcdf,      ONLY: t_input_file
  USE mo_jsb_time,           ONLY: is_newyear, get_year, is_time_ltrig_rad_m1


  ! Use of processes in this module
  dsl4jsb_Use_processes  DYNVEG_, CARBON_, PHENO_, ASSIMI_, DISTURB_, RAD_, TURB_

  ! Use of process configurations
  dsl4jsb_Use_config(DYNVEG_)
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_config(ASSIMI_)
  dsl4jsb_Use_config(DISTURB_)
  dsl4jsb_Use_config(RAD_)
  dsl4jsb_Use_config(TURB_)

  ! Use of process memories
  dsl4jsb_Use_memory(DYNVEG_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(DISTURB_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(TURB_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_dynveg_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dynveg_interface'

  !> Type definition for update_dynveg task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_dynveg
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dynveg    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_dynveg !< Aggregates computed task variables
  END TYPE   tsk_dynveg

  !> Constructor interface for update_dynveg task
  INTERFACE   tsk_dynveg
    PROCEDURE Create_task_dynveg         !< Constructor function for task
  END INTERFACE   tsk_dynveg


CONTAINS
  ! ================================================================================================================================
  !! Constructors for tasks

  !> Constructor for update_dynveg task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "update_dynveg"
  !!
  FUNCTION Create_task_dynveg(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_dynveg::return_ptr)
    CALL return_ptr%Construct(name='dynveg', process_id=DYNVEG_, owner_model_id=model_id)

  END FUNCTION Create_task_dynveg



  ! =======================================================================================================
  !> Register tasks for dynveg process
  !!
  !! @param[in,out] this      Instance of dynveg process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_dynveg_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_dynveg(model_id))
  END SUBROUTINE Register_dynveg_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to update the dynveg state
  !! Task "dynveg" calculates the dynveg state of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE update_dynveg(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER             :: model
    dsl4jsb_Def_config(DYNVEG_)

    ! Local variables
    INTEGER                     :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_dynveg'

    ! Set local variables from options argument
    iblk    = options%iblk

    IF (.NOT. tile%Is_process_active(DYNVEG_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF (.NOT. tile%is_vegetation) THEN
       CALL finish(TRIM(routine), 'Error: JSBACH4 tries to calculate dynveg not on the veg tile!')
    END IF

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(DYNVEG_)

    SELECT CASE (TRIM(dsl4jsb_Config(DYNVEG_)%scheme))
      CASE ('maps')
        CALL dynveg_based_on_maps(tile, options)
!       CASE ('carbon')
!         CALL finish('Scheme carbon/jedi not ready yet')
!         CALL dynveg_based_on_carbon(tile, options)
!       CASE ('jsbach3')
!         CALL finish('Scheme jsbach3 not ready yet')
!         CALL dynveg_based_on_jsbach3(tile, options)
      CASE DEFAULT
        CALL finish('Invalid scheme for dynveg')
    END SELECT

  END SUBROUTINE update_dynveg

  ! ================================================================================================================================
  !>
  !! @brief Implementation of 'MAPS' dynveg scheme
  !!
  !! @param[in,out] tile    Tile for which routine is executed (only for veg-tile)
  !! @param[in]     options Additional run-time parameters.
  !!
  !
  SUBROUTINE dynveg_based_on_maps(tile, options)

    ! Declarations
    USE mo_jsb_time,           ONLY: is_time_experiment_start, is_newyear
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_dynveg_process,     ONLY: relocate_vars_to_new_pft_fractions
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid
    USE mo_jsb_grid,           ONLY: Get_grid
    !USE mo_jsb_memory_class,   ONLY: t_jsb_memory
    !USE mo_jsb_model_usecases,  ONLY: npft

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(DYNVEG_)
    dsl4jsb_Def_memory(DYNVEG_)
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_config(DISTURB_)
    dsl4jsb_Def_memory(DISTURB_)
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_config(TURB_)
    dsl4jsb_Def_memory(TURB_)

    ! Local variables
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    INTEGER                     :: iblk, ics, ice
    INTEGER                     :: current_year
    INTEGER                     :: npft, ipft, iarea, i, j, k, ncanopy, icanopy
    LOGICAL                     :: read_maps = .FALSE.
    TYPE(t_input_file)          :: input_file
    CHARACTER(len=70)           :: filename_maps, varname_within_file
    CHARACTER(len=*), PARAMETER :: routine = modname//':dynveg_based_on_maps'
    CHARACTER(len=70)           :: ipftstring

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & d_fract_pfts,           &
      & d_fract_grow_pfts,      &
      & d_fract_shrink_pfts,    &
      & fract_new_grow_pfts,    &
      & f_old_to_new_grow_pfts, &
      & f_d_to_new_grow_pfts

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: &
      & fract_lost_sum

    ! Declare pointers
    REAL(wp),                    POINTER :: fract_pfts_before_dynvegstep(:,:,:)
    REAL(wp),                    POINTER :: fract_pfts_currentyear_maps(:,:,:)
    !dsl4jsb_Real3D_onChunk :: fract_pfts_before_dynvegstep
    !dsl4jsb_Real3D_onChunk :: fract_pfts_currentyear_maps
    REAL(wp),                    POINTER :: ptr_2D(:,:)
    CLASS(t_jsb_tile_abstract),  POINTER :: ptr_tile

    ! CARBON pfts
    ! Note: all C pool variables and pointers are used for both, _ca and _ta variables!!!
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: &
      & C_pools_sum_before,        &
      & C_pools_sum_after,         &
      & C_pools_sum_check

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & c_green_pfts,              &
      & c_reserve_pfts,            &
      & c_woods_pfts,              &
      & c_crop_harvest_pfts,       &
      !SIZE CLASS 1
      & c_acid_ag1_pfts,           &
      & c_water_ag1_pfts,          &
      & c_ethanol_ag1_pfts,        &
      & c_nonsoluble_ag1_pfts,     &
      & c_acid_bg1_pfts,           &
      & c_water_bg1_pfts,          &
      & c_ethanol_bg1_pfts,        &
      & c_nonsoluble_bg1_pfts,     &
      & c_humus_1_pfts,            &
      !SIZE CLASS 2
      & c_acid_ag2_pfts,           &
      & c_water_ag2_pfts,          &
      & c_ethanol_ag2_pfts,        &
      & c_nonsoluble_ag2_pfts,     &
      & c_acid_bg2_pfts,           &
      & c_water_bg2_pfts,          &
      & c_ethanol_bg2_pfts,        &
      & c_nonsoluble_bg2_pfts,     &
      & c_humus_2_pfts,            &
      & c_bg_sum_pfts,             &
      & c_ag_sum_1_pfts,           &
      & c_litter_green_ag_pfts,    &
      & c_litter_green_bg_pfts,    &
      & c_litter_wood_ag_pfts,     &
      & c_litter_wood_bg_pfts,     &
      !
      & NPP_pot_sum_pfts,          &
      & GPP_sum_pfts,              &
      & alpha_pfts,                &
      & alpha_sum_pfts,            &
      & alpha_yDayMean_pfts,       &
      & soil_respiration_pfts,     &
      & NPP_flux_correction_pfts,  &
      & NPP_pot_yDayMean_pfts,     &
      & NPP_act_yDayMean_pfts,     &
      & excess_NPP_pfts,           &
      & GPP_yDayMean_pfts,         &
      & root_exudates_pfts,        &
      & cflux_c_green_2_herb_pfts,           &
      & cflux_herb_2_littergreen_pfts,       &
      & cflux_herb_2_atm_pfts,               &
      & cflux_dist_greenreserve_2_soil_pfts, &
      & cflux_dist_woods_2_soil_pfts,        &
      & cflux_fire_all_2_atm_pfts,           &
      & current_max_green_pfts,              &
      & veg_carbon_at_max_green_pfts

    ! ASSIMI pfts
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & gross_assimilation_ca_pfts,   &
      & gross_assimilation_pfts,      &
      & dark_respiration_ca_pfts,     &
      & dark_respiration_pfts,        &
      & NPP_pot_rate_ca_pfts,         &
      & canopy_cond_pfts,   &
      & cover_fract_pot_pfts

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
      & assimi_scaling_fact_cl_pfts,   &
      & canopy_cond_cl_pfts,           &
      & gross_assimilation_cl_pfts,    &
      & dark_respiration_cl_pfts,      &
      & carbox_rate_max_cl_pfts,       &
      & e_transport_rate_max_cl_pfts,  &
      & carbox_rate_cl_pfts,           &
      & e_transport_rate_cl_pfts,      &
      & CO2_conc_leaf_cl_pfts

    ! DISTURB pfts
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & cconservation_fire_pfts,  &
      & cconservation_wind_pfts,  &
      & burned_fract_pfts,        &
      & damaged_fract_pfts

    ! PHENO pfts
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & fract_fpc_pfts,                    &
      !& fract_fpc_max_pfts,                &
      !& fract_forest_pfts,                 & ! dont relocate
      & previous_day_NPP_pot_rate_ca_pfts, &
      !& number_of_individuals_pfts,        & ! dont relocate to new pft fractions. Calculated not linear from veg_carbon_at_max_green!
      !& biomass_per_individual_pfts,       & ! same
      !& cconservation_allom_pfts,          &
      & day_NPP_sum_pfts

    ! RAD pfts
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      !& sw_srf_net_pfts,    &
      !& swvis_srf_net_pfts, &
      !& swnir_srf_net_pfts, &
      !& lw_srf_net_pfts,    &
      !& rad_srf_net_pfts,   &
      !& alb_vis_pfts,       & ! only because alb_nir_soil_nolitter
      !& alb_nir_pfts,       & ! only because alb_nir_soil_nolitter
      ! not pfts specific & alb_vis_snow
      ! not pfts specific & alb_nir_snow
      !& alb_vis_soil_pfts,  & ! only because alb_nir_soil_nolitter
      !& alb_nir_soil_pfts,  & ! only because alb_nir_soil_nolitter
      !& alb_vis_lnd_pfts,   &
      !& alb_nir_lnd_pfts,   &
      !& alb_vis_can_pfts,   &
      !& alb_nir_can_pfts,   &
      ! not pft specific & sky_view_fract,       &
      ! not pft specific & sky_view_fract_stem,  &
      & faPAR_pfts,                 &
      & soil_reflectivity_par_pfts    ! because alb_vis_soil
      ! & apar_nacc_pfts

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
      & faPAR_cl_pfts,             &
      & apar_per_lai_cl_pfts
      ! & apar_per_lai_cl_nacc_pfts, &
      ! & apar_per_lai_cl_tavg_pfts

    ! TURB pfts
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & fact_q_air_pfts,           & ! factor that represents snow-water-cover of the surface,
                                     ! which is less a pft property than a result of hydrology
      & fact_qsat_srf_pfts,        & ! see above
      & fact_qsat_trans_srf_pfts     ! see above



    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: &
      & dummy

    ! CARBON avgs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & c_green_lost_avg,           &
      & c_reserve_lost_avg,         &
      & c_woods_lost_avg,           &
      & c_crop_harvest_lost_avg,    &
      !SIZE CLASS 1
      & c_acid_ag1_lost_avg,        &
      & c_water_ag1_lost_avg,       &
      & c_ethanol_ag1_lost_avg,     &
      & c_nonsoluble_ag1_lost_avg,  &
      & c_acid_bg1_lost_avg,        &
      & c_water_bg1_lost_avg,       &
      & c_ethanol_bg1_lost_avg,     &
      & c_nonsoluble_bg1_lost_avg,  &
      & c_humus_1_lost_avg,         &
      !SIZE CLASS 2
      & c_acid_ag2_lost_avg,        &
      & c_water_ag2_lost_avg,       &
      & c_ethanol_ag2_lost_avg,     &
      & c_nonsoluble_ag2_lost_avg,  &
      & c_acid_bg2_lost_avg,        &
      & c_water_bg2_lost_avg,       &
      & c_ethanol_bg2_lost_avg,     &
      & c_nonsoluble_bg2_lost_avg,  &
      & c_humus_2_lost_avg,         &
      & c_bg_sum_lost_avg,          & ! Sum of belowground Yasso pools  [mol(C)^2/m(canopy)]
      & c_ag_sum_1_lost_avg,        & ! Sum of aboveground green Yasso pools  [mol(C)^2/m(canopy)]
      & c_litter_green_ag_lost_avg, &
      & c_litter_green_bg_lost_avg, &
      & c_litter_wood_ag_lost_avg,  &
      & c_litter_wood_bg_lost_avg,  &
      !
      & NPP_pot_sum_lost_avg,          &
      & GPP_sum_lost_avg,              &
      & alpha_lost_avg,                &
      & alpha_sum_lost_avg,            &
      & alpha_yDayMean_lost_avg,       &
      & soil_respiration_lost_avg,     &
      & NPP_flux_correction_lost_avg,  &
      & NPP_pot_yDayMean_lost_avg,     &
      & NPP_act_yDayMean_lost_avg,     &
      & excess_NPP_lost_avg,           &
      & GPP_yDayMean_lost_avg,         &
      & root_exudates_lost_avg,        &
      & cflux_c_green_2_herb_lost_avg,           &
      & cflux_herb_2_littergreen_lost_avg,       &
      & cflux_herb_2_atm_lost_avg,               &
      & cflux_dist_greenreserve_2_soil_lost_avg, &
      & cflux_dist_woods_2_soil_lost_avg,        &
      & cflux_fire_all_2_atm_lost_avg,           &
      & current_max_green_lost_avg,              &
      & veg_carbon_at_max_green_lost_avg

    ! ASSIMI avgs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & gross_assimilation_ca_lost_avg,      &
      & gross_assimilation_lost_avg,         &
      & dark_respiration_ca_lost_avg,        &
      & dark_respiration_lost_avg,           &
      & NPP_pot_rate_ca_lost_avg,            &
      & canopy_cond_lost_avg,      &
      & cover_fract_pot_lost_avg

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
      & assimi_scaling_fact_cl_lost_avg,   &
      & canopy_cond_cl_lost_avg,           &
      & gross_assimilation_cl_lost_avg,    &
      & dark_respiration_cl_lost_avg,      &
      & carbox_rate_max_cl_lost_avg,       &
      & e_transport_rate_max_cl_lost_avg,  &
      & carbox_rate_cl_lost_avg,           &
      & e_transport_rate_cl_lost_avg,      &
      & CO2_conc_leaf_cl_lost_avg

    ! DISTURB avgs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & cconservation_fire_lost_avg,  &
      & cconservation_wind_lost_avg,  &
      & burned_fract_lost_avg,        &
      & damaged_fract_lost_avg

    ! PHENO avgs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & fract_fpc_lost_avg,                    &
      !& fract_fpc_max_lost_avg,                &
      & previous_day_NPP_pot_rate_ca_lost_avg, &
      !& cconservation_allom_lost_avg,          &
      & day_NPP_sum_lost_avg

    ! RAD avgs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      !& sw_srf_net_lost_avg,    &
      !& swvis_srf_net_lost_avg, &
      !& swnir_srf_net_lost_avg, &
      !& lw_srf_net_lost_avg,    &
      !& rad_srf_net_lost_avg,   &
      !& alb_vis_lost_avg,       & ! only because alb_nir_soil_nolitter
      !& alb_nir_lost_avg,       & ! only because alb_nir_soil_nolitter
      ! not pfts specific & alb_vis_snow
      ! not pfts specific & alb_nir_snow
      !& alb_vis_soil_lost_avg, & ! only because alb_nir_soil_nolitter
      !& alb_nir_soil_lost_avg, & ! only because alb_nir_soil_nolitter
      !& alb_vis_lnd_lost_avg,  &
      !& alb_nir_lnd_lost_avg,  &
      !& alb_vis_can_lost_avg,  &
      !& alb_nir_can_lost_avg,  & 
      ! not pft specific & sky_view_fract,       &
      ! not pft specific & sky_view_fract_stem,  &
      & faPAR_lost_avg,                &
      & soil_reflectivity_par_lost_avg   ! because alb_vis_soil
      ! & apar_nacc_lost_avg

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
      & faPAR_cl_lost_avg,             &
      & apar_per_lai_cl_lost_avg
      ! & apar_per_lai_cl_nacc_lost_avg, &
      ! & apar_per_lai_cl_tavg_lost_avg

    ! TURB avgs
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: &
      & fact_q_air_lost_avg,           &
      & fact_qsat_srf_lost_avg,        &
      & fact_qsat_trans_srf_lost_avg



    ! CARBON pointer
    REAL(wp), POINTER, DIMENSION(:,:) ::  &
      & c_green,           &
      & c_reserve,         &
      & c_woods,           &
      & c_crop_harvest,    &
      !SIZE CLASS 1
      & c_acid_ag1,        &
      & c_water_ag1,       &
      & c_ethanol_ag1,     &
      & c_nonsoluble_ag1,  &
      & c_acid_bg1,        &
      & c_water_bg1,       &
      & c_ethanol_bg1,     &
      & c_nonsoluble_bg1,  &
      & c_humus_1,         &
      !SIZE CLASS 2
      & c_acid_ag2,        &
      & c_water_ag2,       &
      & c_ethanol_ag2,     &
      & c_nonsoluble_ag2,  &
      & c_acid_bg2,        &
      & c_water_bg2,       &
      & c_ethanol_bg2,     &
      & c_nonsoluble_bg2,  &
      & c_humus_2,         &
      & c_bg_sum,          &  ! Sum of belowground Yasso pools  [mol(C)^2/m(canopy)]
      & c_ag_sum_1,        &  ! Sum of aboveground green Yasso pools  [mol(C)^2/m(canopy)]
      & c_litter_green_ag, &
      & c_litter_green_bg, &
      & c_litter_wood_ag,  &
      & c_litter_wood_bg,  &
      !
      & NPP_pot_sum,          &
      & GPP_sum,              &
      & alpha,                &
      & alpha_sum,            &
      & alpha_yDayMean,       &
      & soil_respiration,     &
      & NPP_flux_correction,  &
      & NPP_pot_yDayMean,     &
      & NPP_act_yDayMean,     &
      & excess_NPP,           &
      & GPP_yDayMean,         &
      & root_exudates,        &
      & cflux_c_green_2_herb,           &
      & cflux_herb_2_littergreen,       &
      & cflux_herb_2_atm,               &
      & cflux_dist_greenreserve_2_soil, &
      & cflux_dist_woods_2_soil,        &
      & cflux_fire_all_2_atm,           &
      & current_max_green,              &
      & veg_carbon_at_max_green

    ! ASSIMI pointer
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & gross_assimilation_ca,           &
      & gross_assimilation,              &
      & dark_respiration_ca,             &
      & dark_respiration,                &
      & NPP_pot_rate_ca,                 &
      & canopy_cond,           &
      & cover_fract_pot

    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & assimi_scaling_fact_cl,            &
      & canopy_cond_cl,                    &
      & gross_assimilation_cl,             &
      & dark_respiration_cl,               &
      & carbox_rate_max_cl,                &
      & e_transport_rate_max_cl,           &
      & carbox_rate_cl,                    &
      & e_transport_rate_cl,               &
      & CO2_conc_leaf_cl

    ! DISTURB pointer
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & cconservation_fire,  &
      & cconservation_wind,  &
      & burned_fract,        &
      & damaged_fract

    ! PHENO pointer
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & fract_fpc,                    &
      !& fract_fpc_max,                &
      & previous_day_NPP_pot_rate_ca, &
      !& cconservation_allom,          &
      & day_NPP_sum

    ! RAD pointer
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      !& sw_srf_net,    & ! Note: On the leafs, radiation from atmosphere (e.g. swvis_srf_down) and albedos are
                          !       used by the subroutine calc_radiation_surface_net to calculate sw_srf_net.
                          !       Thus rescaling the value of sw_srf_net guarantees the same absolute amount on
                          !       the grid box, but does not keep it consistent with/adapts it to the altered
                          !       surface properties (like albedo)! As we do dynveg directly before the radiation
                          !       time step of the atmos. we want to keep it konsistent with the surface properties!
      !& swvis_srf_net, &
      !& swnir_srf_net, &
      !& lw_srf_net,    &
      !& rad_srf_net,   &
      !& alb_vis,       & ! only because alb_nir_soil_nolitter
      !& alb_nir,       & ! only because alb_nir_soil_nolitter
      ! not pfts specific & alb_vis_snow
      ! not pfts specific & alb_nir_snow
      !& alb_vis_soil, & ! only because alb_nir_soil_nolitter
      !& alb_nir_soil, & ! only because alb_nir_soil_nolitter
      !& alb_vis_lnd,  &
      !& alb_nir_lnd,  &
      !& alb_vis_can,  &
      !& alb_nir_can,  & 
      ! not pft specific & sky_view_fract,       &
      ! not pft specific & sky_view_fract_stem,  &
      & faPAR,                 & ! Note: in contrast to sw_srf_net the par variables shall not be consistent with
                                 ! the surface properties, but shall keep the absolute amount on the grid box
                                 ! because these variables must be consistent with NPP.
      & soil_reflectivity_par    ! because alb_vis_soil
      ! & apar_nacc

    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & faPAR_cl,             &
      & apar_per_lai_cl
      ! & apar_per_lai_cl_nacc, &
      ! & apar_per_lai_cl_tavg

    ! TURB pointer
    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & fact_q_air,         &
      & fact_qsat_srf,      &
      & fact_qsat_trans_srf



    ! Set pointers
    model => Get_model(tile%owner_model_id)
    hgrid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(DYNVEG_)
    dsl4jsb_Get_memory(DYNVEG_)
    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(DISTURB_)
    dsl4jsb_Get_config(RAD_)
    dsl4jsb_Get_config(TURB_)

    fract_pfts_before_dynvegstep => dynveg__mem%fract_pfts_before_dynvegstep ! in
    fract_pfts_currentyear_maps  => dynveg__mem%fract_pfts_currentyear_maps  ! in

    npft    = tile%Get_no_of_children() ! R: Note, we yet have only 11 pfts in the maps input files.
    ncanopy = dsl4jsb_Config(ASSIMI_)%ncanopy

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice


    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(DYNVEG_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! Go

    ! Reading of MAPS for this year
    ! R: Note, we could shift this to SUBROUTINE jsbach_start_timestep, but I want to keep the code "readable".

    ! Read maps once at first day of the year and only in the time step before the radiation time step and only once for the whole grid
    IF (is_newyear(options%current_datetime, options%dtime) .AND. iblk == 1) read_maps = .TRUE.

    IF (read_maps .AND. iblk == 1 .AND. &
      & is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id) ) THEN

      read_maps = .FALSE.

      ! Allocate variables
      ALLOCATE(                                                    &
        & d_fract_pfts(hgrid%nproma, hgrid%nblks, npft),           &
        & d_fract_grow_pfts(hgrid%nproma, hgrid%nblks, npft),      &
        & d_fract_shrink_pfts(hgrid%nproma, hgrid%nblks, npft),    &
        & fract_new_grow_pfts(hgrid%nproma, hgrid%nblks, npft),    &
        & f_old_to_new_grow_pfts(hgrid%nproma, hgrid%nblks, npft), &
        & f_d_to_new_grow_pfts(hgrid%nproma, hgrid%nblks, npft)  )

      ALLOCATE( fract_lost_sum(hgrid%nproma, hgrid%nblks) )
      ALLOCATE( dummy(hgrid%nproma, hgrid%nblks) )

      IF (dsl4jsb_Config(CARBON_)%active == .TRUE.) THEN
        ALLOCATE(                                           &
          & C_pools_sum_before(hgrid%nproma, hgrid%nblks),  &
          & C_pools_sum_after(hgrid%nproma, hgrid%nblks),   &
          ! & C_pools_sum_check(hgrid%nproma, hgrid%nblks),   &
          !
          !
          & c_green_pfts(hgrid%nproma, hgrid%nblks, npft),              &
          & c_reserve_pfts(hgrid%nproma, hgrid%nblks, npft),            &
          & c_woods_pfts(hgrid%nproma, hgrid%nblks, npft),              &
          & c_crop_harvest_pfts(hgrid%nproma, hgrid%nblks, npft),       &
          !SIZE CLASS 1
          & c_acid_ag1_pfts(hgrid%nproma, hgrid%nblks, npft),           &
          & c_water_ag1_pfts(hgrid%nproma, hgrid%nblks, npft),          &
          & c_ethanol_ag1_pfts(hgrid%nproma, hgrid%nblks, npft),        &
          & c_nonsoluble_ag1_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & c_acid_bg1_pfts(hgrid%nproma, hgrid%nblks, npft),           &
          & c_water_bg1_pfts(hgrid%nproma, hgrid%nblks, npft),          &
          & c_ethanol_bg1_pfts(hgrid%nproma, hgrid%nblks, npft),        &
          & c_nonsoluble_bg1_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & c_humus_1_pfts(hgrid%nproma, hgrid%nblks, npft),            &
          !SIZE CLASS 2
          & c_acid_ag2_pfts(hgrid%nproma, hgrid%nblks, npft),           &
          & c_water_ag2_pfts(hgrid%nproma, hgrid%nblks, npft),          &
          & c_ethanol_ag2_pfts(hgrid%nproma, hgrid%nblks, npft),        &
          & c_nonsoluble_ag2_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & c_acid_bg2_pfts(hgrid%nproma, hgrid%nblks, npft),           &
          & c_water_bg2_pfts(hgrid%nproma, hgrid%nblks, npft),          &
          & c_ethanol_bg2_pfts(hgrid%nproma, hgrid%nblks, npft),        &
          & c_nonsoluble_bg2_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & c_humus_2_pfts(hgrid%nproma, hgrid%nblks, npft),            &
          !
          & c_bg_sum_pfts(hgrid%nproma, hgrid%nblks, npft),             &   ! Sum of belowground Yasso pools  [mol(C)^2/m(canopy)]
          & c_ag_sum_1_pfts(hgrid%nproma, hgrid%nblks, npft),           &   ! Sum of aboveground green Yasso pools  [mol(C)^2/m(canopy)]
          & c_litter_green_ag_pfts(hgrid%nproma, hgrid%nblks, npft),    &
          & c_litter_green_bg_pfts(hgrid%nproma, hgrid%nblks, npft),    &
          & c_litter_wood_ag_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & c_litter_wood_bg_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          !
          & NPP_pot_sum_pfts(hgrid%nproma, hgrid%nblks, npft),          &
          & GPP_sum_pfts(hgrid%nproma, hgrid%nblks, npft),              &
          & alpha_pfts(hgrid%nproma, hgrid%nblks, npft),                &
          & alpha_sum_pfts(hgrid%nproma, hgrid%nblks, npft),            &
          & alpha_yDayMean_pfts(hgrid%nproma, hgrid%nblks, npft),       &
          & soil_respiration_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & NPP_flux_correction_pfts(hgrid%nproma, hgrid%nblks, npft),  &
          & NPP_pot_yDayMean_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & NPP_act_yDayMean_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & excess_NPP_pfts(hgrid%nproma, hgrid%nblks, npft),           &
          & GPP_yDayMean_pfts(hgrid%nproma, hgrid%nblks, npft),         &
          & root_exudates_pfts(hgrid%nproma, hgrid%nblks, npft),        &
          & cflux_c_green_2_herb_pfts(hgrid%nproma, hgrid%nblks, npft),           &
          & cflux_herb_2_littergreen_pfts(hgrid%nproma, hgrid%nblks, npft),       &
          & cflux_herb_2_atm_pfts(hgrid%nproma, hgrid%nblks, npft),               &
          & cflux_dist_greenreserve_2_soil_pfts(hgrid%nproma, hgrid%nblks, npft), &
          & cflux_dist_woods_2_soil_pfts(hgrid%nproma, hgrid%nblks, npft),        &
          & cflux_fire_all_2_atm_pfts(hgrid%nproma, hgrid%nblks, npft),           &
          & current_max_green_pfts(hgrid%nproma, hgrid%nblks, npft),              &
          & veg_carbon_at_max_green_pfts(hgrid%nproma, hgrid%nblks, npft)         &
            )


        ALLOCATE(                                                        &
          & c_green_lost_avg(hgrid%nproma, hgrid%nblks, npft),           &
          & c_reserve_lost_avg(hgrid%nproma, hgrid%nblks, npft),         &
          & c_woods_lost_avg(hgrid%nproma, hgrid%nblks, npft),           &
          & c_crop_harvest_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          !SIZE CLASS 1
          & c_acid_ag1_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &
          & c_water_ag1_lost_avg(hgrid%nproma, hgrid%nblks, npft),       &
          & c_ethanol_ag1_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          & c_nonsoluble_ag1_lost_avg(hgrid%nproma, hgrid%nblks, npft),  &
          & c_acid_bg1_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &
          & c_water_bg1_lost_avg(hgrid%nproma, hgrid%nblks, npft),       &
          & c_ethanol_bg1_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          & c_nonsoluble_bg1_lost_avg(hgrid%nproma, hgrid%nblks, npft),  &
          & c_humus_1_lost_avg(hgrid%nproma, hgrid%nblks, npft),         &
          !SIZE CLASS 2
          & c_acid_ag2_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &
          & c_water_ag2_lost_avg(hgrid%nproma, hgrid%nblks, npft),       &
          & c_ethanol_ag2_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          & c_nonsoluble_ag2_lost_avg(hgrid%nproma, hgrid%nblks, npft),  &
          & c_acid_bg2_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &
          & c_water_bg2_lost_avg(hgrid%nproma, hgrid%nblks, npft),       &
          & c_ethanol_bg2_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          & c_nonsoluble_bg2_lost_avg(hgrid%nproma, hgrid%nblks, npft),  &
          & c_humus_2_lost_avg(hgrid%nproma, hgrid%nblks, npft),         &
          & c_bg_sum_lost_avg(hgrid%nproma, hgrid%nblks, npft),          &  ! Sum of belowground Yasso pools  [mol(C)^2/m(canopy)]
          & c_ag_sum_1_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &  ! Sum of aboveground green Yasso pools  [mol(C)^2/m(canopy)]
          & c_litter_green_ag_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          & c_litter_green_bg_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          & c_litter_wood_ag_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          & c_litter_wood_bg_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          !
          & NPP_pot_sum_lost_avg(hgrid%nproma, hgrid%nblks, npft),          &
          & GPP_sum_lost_avg(hgrid%nproma, hgrid%nblks, npft),              &
          & alpha_lost_avg(hgrid%nproma, hgrid%nblks, npft),                &
          & alpha_sum_lost_avg(hgrid%nproma, hgrid%nblks, npft),            &
          & alpha_yDayMean_lost_avg(hgrid%nproma, hgrid%nblks, npft),       &
          & soil_respiration_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          & NPP_flux_correction_lost_avg(hgrid%nproma, hgrid%nblks, npft),  &
          & NPP_pot_yDayMean_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          & NPP_act_yDayMean_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          & excess_NPP_lost_avg(hgrid%nproma, hgrid%nblks, npft),           &
          & GPP_yDayMean_lost_avg(hgrid%nproma, hgrid%nblks, npft),         &
          & root_exudates_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &
          & cflux_c_green_2_herb_lost_avg(hgrid%nproma, hgrid%nblks, npft),           &
          & cflux_herb_2_littergreen_lost_avg(hgrid%nproma, hgrid%nblks, npft),       &
          & cflux_herb_2_atm_lost_avg(hgrid%nproma, hgrid%nblks, npft),               &
          & cflux_dist_greenreserve_2_soil_lost_avg(hgrid%nproma, hgrid%nblks, npft), &
          & cflux_dist_woods_2_soil_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &
          & cflux_fire_all_2_atm_lost_avg(hgrid%nproma, hgrid%nblks, npft),           &
          & current_max_green_lost_avg(hgrid%nproma, hgrid%nblks, npft),              &
          & veg_carbon_at_max_green_lost_avg(hgrid%nproma, hgrid%nblks, npft)         &
          )
      END IF ! CARBON is active


      IF (dsl4jsb_Config(ASSIMI_)%active == .TRUE.) THEN
        ALLOCATE(                                                          &
          & gross_assimilation_ca_pfts(hgrid%nproma, hgrid%nblks, npft),   &
          & gross_assimilation_pfts(hgrid%nproma, hgrid%nblks, npft),      &
          & dark_respiration_ca_pfts(hgrid%nproma, hgrid%nblks, npft),     &
          & dark_respiration_pfts(hgrid%nproma, hgrid%nblks, npft),        &
          & NPP_pot_rate_ca_pfts(hgrid%nproma, hgrid%nblks, npft),         &
          & canopy_cond_pfts(hgrid%nproma, hgrid%nblks, npft),   &
          & cover_fract_pot_pfts(hgrid%nproma, hgrid%nblks, npft)          &
          & )

        ALLOCATE(                                                                    &
          & assimi_scaling_fact_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),   &
          & canopy_cond_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),           &
          & gross_assimilation_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),    &
          & dark_respiration_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),      &
          & carbox_rate_max_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),       &
          & e_transport_rate_max_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),  &
          & carbox_rate_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),           &
          & e_transport_rate_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),      &
          & CO2_conc_leaf_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft)          &
          & )

        ALLOCATE(                                                            &
          & gross_assimilation_ca_lost_avg(hgrid%nproma, hgrid%nblks, npft), &
          & gross_assimilation_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          & dark_respiration_ca_lost_avg(hgrid%nproma, hgrid%nblks, npft),   &
          & dark_respiration_lost_avg(hgrid%nproma, hgrid%nblks, npft),      &
          & NPP_pot_rate_ca_lost_avg(hgrid%nproma, hgrid%nblks, npft),       &
          & canopy_cond_lost_avg(hgrid%nproma, hgrid%nblks, npft), &
          & cover_fract_pot_lost_avg(hgrid%nproma, hgrid%nblks, npft)        &
          & )

        ALLOCATE(                                                                        &
          & assimi_scaling_fact_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),   &
          & canopy_cond_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),           &
          & gross_assimilation_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),    &
          & dark_respiration_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),      &
          & carbox_rate_max_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),       &
          & e_transport_rate_max_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),  &
          & carbox_rate_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),           &
          & e_transport_rate_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),      &
          & CO2_conc_leaf_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft)          &
          & )
      END IF ! ASSIMI is active


      IF (dsl4jsb_Config(DISTURB_)%active == .TRUE.) THEN
        ALLOCATE(                                                     &
          & cconservation_fire_pfts(hgrid%nproma, hgrid%nblks, npft), &
          & cconservation_wind_pfts(hgrid%nproma, hgrid%nblks, npft), &
          & burned_fract_pfts(hgrid%nproma, hgrid%nblks, npft),       &
          & damaged_fract_pfts(hgrid%nproma, hgrid%nblks, npft)       &
          & )
  
        ALLOCATE(                                                         &
          & cconservation_fire_lost_avg(hgrid%nproma, hgrid%nblks, npft), &
          & cconservation_wind_lost_avg(hgrid%nproma, hgrid%nblks, npft), &
          & burned_fract_lost_avg(hgrid%nproma, hgrid%nblks, npft),       &
          & damaged_fract_lost_avg(hgrid%nproma, hgrid%nblks, npft)       &
          & )
      END IF ! DISTURB is active


      IF (dsl4jsb_Config(PHENO_)%active == .TRUE.) THEN
        ALLOCATE(                                                               &
          & fract_fpc_pfts(hgrid%nproma, hgrid%nblks, npft),                    &
          !& fract_fpc_max_pfts(hgrid%nproma, hgrid%nblks, npft),                &
          & previous_day_NPP_pot_rate_ca_pfts(hgrid%nproma, hgrid%nblks, npft), &
          !IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
          !& cconservation_allom_pfts(hgrid%nproma, hgrid%nblks, npft),          &
          & day_NPP_sum_pfts(hgrid%nproma, hgrid%nblks, npft)                   &
          & )

        ALLOCATE(                                                                   &
          & fract_fpc_lost_avg(hgrid%nproma, hgrid%nblks, npft),                    &
          !& fract_fpc_max_lost_avg(hgrid%nproma, hgrid%nblks, npft),                &
          & previous_day_NPP_pot_rate_ca_lost_avg(hgrid%nproma, hgrid%nblks, npft), &
          !IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
          !& cconservation_allom_lost_avg(hgrid%nproma, hgrid%nblks, npft),          &
          & day_NPP_sum_lost_avg(hgrid%nproma, hgrid%nblks, npft)                   &
          & )
      END IF ! PHENO is active

      IF (dsl4jsb_Config(RAD_)%active == .TRUE.) THEN
        ALLOCATE(                                                &
          !& sw_srf_net_pfts(hgrid%nproma, hgrid%nblks, npft),    &
          !& swvis_srf_net_pfts(hgrid%nproma, hgrid%nblks, npft), &
          !& swnir_srf_net_pfts(hgrid%nproma, hgrid%nblks, npft), &
          !& lw_srf_net_pfts(hgrid%nproma, hgrid%nblks, npft),    &
          !& rad_srf_net_pfts(hgrid%nproma, hgrid%nblks, npft),   &
          !& alb_vis_pfts(hgrid%nproma, hgrid%nblks, npft),       &
          !& alb_nir_pfts(hgrid%nproma, hgrid%nblks, npft),       &
          !& alb_vis_soil_pfts(hgrid%nproma, hgrid%nblks, npft),  &
          !& alb_nir_soil_pfts(hgrid%nproma, hgrid%nblks, npft),  &
          !& alb_vis_lnd_pfts(hgrid%nproma, hgrid%nblks, npft),   &
          !& alb_nir_lnd_pfts(hgrid%nproma, hgrid%nblks, npft),   &
          !& alb_vis_can_pfts(hgrid%nproma, hgrid%nblks, npft),   &
          !& alb_nir_can_pfts(hgrid%nproma, hgrid%nblks, npft),   &
          & faPAR_pfts(hgrid%nproma, hgrid%nblks, npft),                 &
          & soil_reflectivity_par_pfts(hgrid%nproma, hgrid%nblks, npft)  &
          ! & apar_nacc_pfts(hgrid%nproma, hgrid%nblks, npft)              &
          & )

        ALLOCATE(                                                                &
          & faPAR_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft),             &
          & apar_per_lai_cl_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft)       &
          ! & apar_per_lai_cl_nacc_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft), &
          ! & apar_per_lai_cl_tavg_pfts(hgrid%nproma, ncanopy, hgrid%nblks, npft)  &
          & )

        ALLOCATE(                                                     &
          !& sw_srf_net_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          !& swvis_srf_net_lost_avg(hgrid%nproma, hgrid%nblks, npft),  &
          !& swnir_srf_net_lost_avg(hgrid%nproma, hgrid%nblks, npft),  &
          !& lw_srf_net_lost_avg(hgrid%nproma, hgrid%nblks, npft),     &
          !& rad_srf_net_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          !& alb_vis_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &
          !& alb_nir_lost_avg(hgrid%nproma, hgrid%nblks, npft),        &
          !& alb_vis_soil_lost_avg(hgrid%nproma, hgrid%nblks, npft),   &
          !& alb_nir_soil_lost_avg(hgrid%nproma, hgrid%nblks, npft),   &
          !& alb_vis_lnd_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          !& alb_nir_lnd_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          !& alb_vis_can_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          !& alb_nir_can_lost_avg(hgrid%nproma, hgrid%nblks, npft),    &
          & faPAR_lost_avg(hgrid%nproma, hgrid%nblks, npft),                 &
          & soil_reflectivity_par_lost_avg(hgrid%nproma, hgrid%nblks, npft)  &
          ! & apar_nacc_lost_avg(hgrid%nproma, hgrid%nblks, npft)              &
          & )

        ALLOCATE(                                                                    &
          & faPAR_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft),             &
          & apar_per_lai_cl_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft)       &
          ! & apar_per_lai_cl_nacc_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft), &
          ! & apar_per_lai_cl_tavg_lost_avg(hgrid%nproma, ncanopy, hgrid%nblks, npft)  &
          & )
      END IF ! RAD is active


      IF (dsl4jsb_Config(TURB_)%active == .TRUE.) THEN
        ALLOCATE(                                                       &
          & fact_q_air_pfts(hgrid%nproma, hgrid%nblks, npft),           &
          & fact_qsat_srf_pfts(hgrid%nproma, hgrid%nblks, npft),        &
          & fact_qsat_trans_srf_pfts(hgrid%nproma, hgrid%nblks, npft)   & 
          & )

        ALLOCATE(                                                         &
          & fact_q_air_lost_avg(hgrid%nproma, hgrid%nblks, npft),         &
          & fact_qsat_srf_lost_avg(hgrid%nproma, hgrid%nblks, npft),      &
          & fact_qsat_trans_srf_lost_avg(hgrid%nproma, hgrid%nblks, npft) &
          & )
      END IF ! TURB is active


      ! Generate name of map file with new landcover, e.g. cover_fract_T63GR15_11tiles_1991.nc
      current_year  = INT(get_year(options%current_datetime))
      IF (current_year < 10000) THEN
        ! The filename is fixed in JS3. In JS4 we could read it from namelist, if worth it.
        WRITE(filename_maps,'(a,I4.4,a)') "bc_land_frac_11pfts_", current_year, ".nc"
      ELSE
        CALL finish('read_landcover_fractions from map file','Only years between 0 and 9999 supported currently')
      END IF

      input_file = jsb_netcdf_open_input(TRIM(filename_maps), model%grid_id)

      DO ipft=1,npft

        ! Point to the correct pft
        IF (ipft == 1) THEN
          ptr_tile => tile%Get_first_child_tile()
        ELSE
          ptr_tile => ptr_tile%Get_next_sibling_tile()
        END IF

        ! If necessary collect the actual pft-fractions to remember the fractions before this dynveg step
        IF (is_time_experiment_start(options%current_datetime)) THEN  ! Start of experiment
            !   We yet do not have reasonable values for fract_pfts_before_dynvegstep
            !   (except fract_pfts_before_dynvegstep was included in a restart file which is used for this experiment).
            !   Therefore we take the actual fractions (received from initialization or restart file)
            CALL ptr_tile%Get_fraction(fract=fract_pfts_before_dynvegstep(:,:,ipft))
        !ELSE we have to do nothing because we just take the existing fract_pfts_before_dynvegstep
        END IF

        ! Read this years pft-fraction from map file
        write(ipftstring , *) ipft
        ! IF (debug_on()) CALL message(TRIM(routine), '...reading yearly cover_fract maps on pft'//TRIM(ipftstring)//' ...')
        CALL message(TRIM(routine), '...reading yearly cover_fract maps on pft'//TRIM(ipftstring)//' ...')
        WRITE(varname_within_file,'(a,I2.2)') "fract_pft", ipft

        ptr_2D => input_file%Read_2d(             &
          & variable_name=varname_within_file,    &
          & fill_array = dynveg__mem%fract_pfts_currentyear_maps(:,:,ipft) )

        ! Avoid negative values and inconsistencies  
        dynveg__mem%fract_pfts_currentyear_maps(:,:,ipft) =                    &
          & MERGE(dynveg__mem%fract_pfts_currentyear_maps(:,:,ipft),           &
          & 0._wp, dynveg__mem%fract_pfts_currentyear_maps(:,:,ipft) >= 0._wp)

          
          

        ! Set pft fraction according to the map
        CALL ptr_tile%Set_fraction(fract=dynveg__mem%fract_pfts_currentyear_maps(:,:,ipft))
        ! Note, all other tile fractions should not be affected by the pft fraction changes!
        
      END DO ! ipft

      ! As soon as we will get absolute instead of relative fractions in JS4 we have to aggregate the pft fraction changes!
      ! However, at that time we have to consider, that the pfts in the maps-file sum up to 1, which does only agree with
      ! our usecase JS4 tileing if we interpret the map-fractions as relative! In this case we have to convert them to absolute
      ! values!
      ! As the pfts in the maps-file are regarded relative and the veg tile and the veg tile fraction of 1 remains constant
      ! for each year, we do not want to set the veg tile here!
      ! => we neither aggregate the pft fractions nor set the veg tile fractions


      ! Recalculate those variables that differ between pfts

      ! Calculate needed delta fractions
      d_fract_pfts(:,:,:)        = fract_pfts_currentyear_maps(:,:,:) - fract_pfts_before_dynvegstep(:,:,:)
      d_fract_grow_pfts(:,:,:)   = MAX(0.0_wp, d_fract_pfts(:,:,:))             ! select deltas of pfts which grow (for each gridcell)
      d_fract_shrink_pfts(:,:,:) = MAX(0.0_wp, d_fract_pfts(:,:,:) * (-1._wp))  ! same for the shinking pfts
                                                                                ! MIN(0.0_wp, d_fract_pfts(:,:)) would be negative!
      ! Calculate often used fraction factors
      fract_new_grow_pfts(:,:,:)    = fract_pfts_before_dynvegstep(:,:,:) + d_fract_grow_pfts(:,:,:)
      fract_lost_sum(:,:)           = SUM(d_fract_shrink_pfts(:,:,:),DIM=3)
      f_old_to_new_grow_pfts(:,:,:) = fract_pfts_before_dynvegstep(:,:,:) / fract_new_grow_pfts(:,:,:)
      f_d_to_new_grow_pfts(:,:,:)   = d_fract_grow_pfts(:,:,:)            / fract_new_grow_pfts(:,:,:)


      ! We have to handle divisions by zero in a way that avoids special treatments in the code:
      ! (1) Divisions by zero would happen when there are no fraction changes on all pfts of the veg tile.
      !     On these grid points we want the *_lost_avg variables to be small to guarantee no "real" changes
      WHERE ( fract_lost_sum(:,:) < EPSILON(1._wp) )
          fract_lost_sum(:,:) = 1.e16_wp ! We just use a arbitrarily high value to guarantee *_lost_avg gets small
      END WHERE
      ! (2) Divisions by zero happened when fract_new_grow_pfts(:,:,:) where zero.
      !     Therefore we have to correct the resulting NaN values in a way to keep the old variable values (no "real" changes)
      WHERE ( fract_new_grow_pfts(:,:,:) <EPSILON(1._wp) )
          f_old_to_new_grow_pfts(:,:,:) = 1._wp
          f_d_to_new_grow_pfts(:,:,:)   = 0._wp
      END WHERE


      IF (dsl4jsb_Config(CARBON_)%active == .TRUE.) THEN
        IF (debug_on()) CALL message(TRIM(routine), '...starting to recalculate carbon variables according to '//TRIM(filename_maps)//' ...')

        DO iarea=1,2 ! because we want to do the same for _ca and _ta variables!!!
          ! Collect pfts of the variables
          ! R: mo_jsb_memory_class.f90: FUNCTION Collect_var_real3d could also be used.
          !    Not used because ics, ice, iblk are used there.
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(CARBON_)%p)
                  CLASS IS (t_carbon_memory)
                    carbon__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_carbon_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            IF (iarea == 1 ) THEN ! _ca variables
              c_green           => carbon__mem%c_green%ptr(:,:)           ! inout
              c_reserve         => carbon__mem%c_reserve%ptr(:,:)         ! inout
              c_woods           => carbon__mem%c_woods%ptr(:,:)           ! inout
              c_crop_harvest    => carbon__mem%c_crop_harvest%ptr(:,:)    ! inout
              c_litter_green_ag => carbon__mem%c_litter_green_ag%ptr(:,:) ! inout
              c_litter_green_bg => carbon__mem%c_litter_green_bg%ptr(:,:) ! inout
              c_litter_wood_ag  => carbon__mem%c_litter_wood_ag %ptr(:,:) ! inout
              c_litter_wood_bg  => carbon__mem%c_litter_wood_bg %ptr(:,:) ! inout
              ! size class 1
              c_acid_ag1        => carbon__mem%c_acid_ag1%ptr(:,:)        ! inout
              c_water_ag1       => carbon__mem%c_water_ag1%ptr(:,:)       ! inout
              c_ethanol_ag1     => carbon__mem%c_ethanol_ag1%ptr(:,:)     ! inout
              c_nonsoluble_ag1  => carbon__mem%c_nonsoluble_ag1%ptr(:,:)  ! inout
              c_acid_bg1        => carbon__mem%c_acid_bg1%ptr(:,:)        ! inout
              c_water_bg1       => carbon__mem%c_water_bg1%ptr(:,:)       ! inout
              c_ethanol_bg1     => carbon__mem%c_ethanol_bg1%ptr(:,:)     ! inout
              c_nonsoluble_bg1  => carbon__mem%c_nonsoluble_bg1%ptr(:,:)  ! inout
              c_humus_1         => carbon__mem%c_humus_1%ptr(:,:)         ! inout
              ! size class 2
              c_acid_ag2        => carbon__mem%c_acid_ag2%ptr(:,:)        ! inout
              c_water_ag2       => carbon__mem%c_water_ag2%ptr(:,:)       ! inout
              c_ethanol_ag2     => carbon__mem%c_ethanol_ag2%ptr(:,:)     ! inout
              c_nonsoluble_ag2  => carbon__mem%c_nonsoluble_ag2%ptr(:,:)  ! inout
              c_acid_bg2        => carbon__mem%c_acid_bg2%ptr(:,:)        ! inout
              c_water_bg2       => carbon__mem%c_water_bg2%ptr(:,:)       ! inout
              c_ethanol_bg2     => carbon__mem%c_ethanol_bg2%ptr(:,:)     ! inout
              c_nonsoluble_bg2  => carbon__mem%c_nonsoluble_bg2%ptr(:,:)  ! inout
              c_humus_2         => carbon__mem%c_humus_2%ptr(:,:)         ! inout
              c_bg_sum          => carbon__mem%c_bg_sum%ptr(:,:)          ! inout
              c_ag_sum_1        => carbon__mem%c_ag_sum_1%ptr(:,:)        ! inout
              c_litter_green_ag => carbon__mem%c_litter_green_ag%ptr(:,:) ! inout
              c_litter_green_bg => carbon__mem%c_litter_green_bg%ptr(:,:) ! inout
              c_litter_wood_ag  => carbon__mem%c_litter_wood_ag%ptr(:,:)  ! inout
              c_litter_wood_bg  => carbon__mem%c_litter_wood_bg%ptr(:,:)  ! inout
              !
              NPP_pot_sum          => carbon__mem%NPP_pot_sum%ptr(:,:)    ! inout
              GPP_sum              => carbon__mem%GPP_sum%ptr(:,:)        ! inout
              alpha                => carbon__mem%alpha%ptr(:,:)          ! inout
              alpha_sum            => carbon__mem%alpha_sum%ptr(:,:)      ! inout
              alpha_yDayMean       => carbon__mem%alpha_yDayMean%ptr(:,:) ! inout
              soil_respiration     => carbon__mem%soil_respiration%ptr(:,:)    ! inout
              NPP_flux_correction  => carbon__mem%NPP_flux_correction%ptr(:,:) ! inout
              NPP_pot_yDayMean     => carbon__mem%NPP_pot_yDayMean%ptr(:,:)    ! inout
              NPP_act_yDayMean     => carbon__mem%NPP_act_yDayMean%ptr(:,:)    ! inout
              excess_NPP           => carbon__mem%excess_NPP%ptr(:,:)          ! inout
              GPP_yDayMean         => carbon__mem%GPP_yDayMean%ptr(:,:)        ! inout
              root_exudates        => carbon__mem%root_exudates%ptr(:,:)       ! inout
              cflux_c_green_2_herb => carbon__mem%cflux_c_green_2_herb%ptr(:,:)         ! inout
              cflux_herb_2_littergreen => carbon__mem%cflux_herb_2_littergreen%ptr(:,:) ! inout
              cflux_herb_2_atm     => carbon__mem%cflux_herb_2_atm%ptr(:,:)             ! inout
              cflux_dist_greenreserve_2_soil => carbon__mem%cflux_dist_greenreserve_2_soil%ptr(:,:) ! inout
              cflux_dist_woods_2_soil => carbon__mem%cflux_dist_woods_2_soil%ptr(:,:)   ! inout
              cflux_fire_all_2_atm => carbon__mem%cflux_fire_all_2_atm%ptr(:,:)         ! inout
              IF ( dsl4jsb_Config(PHENO_)%active == .TRUE.  .AND.   dsl4jsb_Config(PHENO_)%l_forestRegrowth ) THEN
                current_max_green       => carbon__mem%current_max_green%ptr(:,:)       ! inout
                veg_carbon_at_max_green => carbon__mem%veg_carbon_at_max_green%ptr(:,:) ! inout
              END IF
            ELSE ! _ta variables
                 ! * Note: the _ta var does not exist, therefore we point it again to the ca var.
              c_green           => carbon__mem%c_green_ta%ptr(:,:)           ! inout
              c_reserve         => carbon__mem%c_reserve_ta%ptr(:,:)         ! inout
              c_woods           => carbon__mem%c_woods_ta%ptr(:,:)           ! inout
              c_crop_harvest    => carbon__mem%c_crop_harvest_ta%ptr(:,:)    ! inout
              c_litter_green_ag => carbon__mem%c_litter_green_ag_ta%ptr(:,:) ! inout
              c_litter_green_bg => carbon__mem%c_litter_green_bg_ta%ptr(:,:) ! inout
              c_litter_wood_ag  => carbon__mem%c_litter_wood_ag_ta%ptr(:,:)  ! inout
              c_litter_wood_bg  => carbon__mem%c_litter_wood_bg_ta%ptr(:,:)  ! inout
              ! size class 1
              c_acid_ag1        => carbon__mem%c_acid_ag1_ta%ptr(:,:)        ! inout
              c_water_ag1       => carbon__mem%c_water_ag1_ta%ptr(:,:)       ! inout
              c_ethanol_ag1     => carbon__mem%c_ethanol_ag1_ta%ptr(:,:)     ! inout
              c_nonsoluble_ag1  => carbon__mem%c_nonsoluble_ag1_ta%ptr(:,:)  ! inout
              c_acid_bg1        => carbon__mem%c_acid_bg1_ta%ptr(:,:)        ! inout
              c_water_bg1       => carbon__mem%c_water_bg1_ta%ptr(:,:)       ! inout
              c_ethanol_bg1     => carbon__mem%c_ethanol_bg1_ta%ptr(:,:)     ! inout
              c_nonsoluble_bg1  => carbon__mem%c_nonsoluble_bg1_ta%ptr(:,:)  ! inout
              c_humus_1         => carbon__mem%c_humus_1_ta%ptr(:,:)         ! inout
              ! size class 2
              c_acid_ag2        => carbon__mem%c_acid_ag2_ta%ptr(:,:)        ! inout
              c_water_ag2       => carbon__mem%c_water_ag2_ta%ptr(:,:)       ! inout
              c_ethanol_ag2     => carbon__mem%c_ethanol_ag2_ta%ptr(:,:)     ! inout
              c_nonsoluble_ag2  => carbon__mem%c_nonsoluble_ag2_ta%ptr(:,:)  ! inout
              c_acid_bg2        => carbon__mem%c_acid_bg2_ta%ptr(:,:)        ! inout
              c_water_bg2       => carbon__mem%c_water_bg2_ta%ptr(:,:)       ! inout
              c_ethanol_bg2     => carbon__mem%c_ethanol_bg2_ta%ptr(:,:)     ! inout
              c_nonsoluble_bg2  => carbon__mem%c_nonsoluble_bg2_ta%ptr(:,:)  ! inout
              c_humus_2         => carbon__mem%c_humus_2_ta%ptr(:,:)         ! inout
              c_bg_sum          => carbon__mem%c_bg_sum%ptr(:,:)             ! inout *
              c_ag_sum_1        => carbon__mem%c_ag_sum_1%ptr(:,:)           ! inout *
              c_litter_green_ag => carbon__mem%c_litter_green_ag_ta%ptr(:,:) ! inout
              c_litter_green_bg => carbon__mem%c_litter_green_bg_ta%ptr(:,:) ! inout
              c_litter_wood_ag  => carbon__mem%c_litter_wood_ag_ta%ptr(:,:)  ! inout
              c_litter_wood_bg  => carbon__mem%c_litter_wood_bg_ta%ptr(:,:)  ! inout
              !
              NPP_pot_sum          => carbon__mem%NPP_pot_sum%ptr(:,:)       ! inout
              GPP_sum              => carbon__mem%GPP_sum%ptr(:,:)           ! inout *
              alpha                => carbon__mem%alpha%ptr(:,:)             ! inout *
              alpha_sum            => carbon__mem%alpha_sum%ptr(:,:)         ! inout *
              alpha_yDayMean       => carbon__mem%alpha_yDayMean%ptr(:,:)    ! inout *
              soil_respiration     => carbon__mem%soil_respiration_ta%ptr(:,:)     ! inout
              NPP_flux_correction  => carbon__mem%NPP_flux_correction_ta%ptr(:,:)  ! inout
              NPP_pot_yDayMean     => carbon__mem%NPP_pot_yDayMean_ta%ptr(:,:)     ! inout
              NPP_act_yDayMean     => carbon__mem%NPP_act_yDayMean_ta%ptr(:,:)     ! inout
              excess_NPP           => carbon__mem%excess_NPP_ta%ptr(:,:)           ! inout
              GPP_yDayMean         => carbon__mem%GPP_yDayMean_ta%ptr(:,:)         ! inout
              root_exudates        => carbon__mem%root_exudates_ta%ptr(:,:)        ! inout
              cflux_c_green_2_herb => carbon__mem%cflux_c_green_2_herb_ta%ptr(:,:) ! inout
              cflux_herb_2_littergreen => carbon__mem%cflux_herb_2_littergreen%ptr(:,:)                ! inout *
              cflux_herb_2_atm     => carbon__mem%cflux_herb_2_atm%ptr(:,:)                            ! inout *
              cflux_dist_greenreserve_2_soil => carbon__mem%cflux_dist_greenreserve_2_soil_ta%ptr(:,:) ! inout
              cflux_dist_woods_2_soil => carbon__mem%cflux_dist_woods_2_soil_ta%ptr(:,:) ! inout
              cflux_fire_all_2_atm => carbon__mem%cflux_fire_all_2_atm_ta%ptr(:,:)       ! inout
              IF ( dsl4jsb_Config(PHENO_)%active == .TRUE.  .AND.   dsl4jsb_Config(PHENO_)%l_forestRegrowth ) THEN
                current_max_green       => carbon__mem%current_max_green%ptr(:,:)          ! inout *
                veg_carbon_at_max_green => carbon__mem%veg_carbon_at_max_green%ptr(:,:)    ! inout *
              END IF
            END IF

            ! Collect
            c_green_pfts(:,:,ipft)           = c_green(:,:)
            c_woods_pfts(:,:,ipft)           = c_woods(:,:)
            c_reserve_pfts(:,:,ipft)         = c_reserve(:,:)
            c_crop_harvest_pfts(:,:,ipft)    = c_crop_harvest(:,:)
            ! size class 1
            c_acid_ag1_pfts(:,:,ipft)        = c_acid_ag1(:,:)
            c_water_ag1_pfts(:,:,ipft)       = c_water_ag1(:,:)
            c_ethanol_ag1_pfts(:,:,ipft)     = c_ethanol_ag1(:,:)
            c_nonsoluble_ag1_pfts(:,:,ipft)  = c_nonsoluble_ag1(:,:)
            c_acid_bg1_pfts(:,:,ipft)        = c_acid_bg1(:,:)
            c_water_bg1_pfts(:,:,ipft)       = c_water_bg1(:,:)
            c_ethanol_bg1_pfts(:,:,ipft)     = c_ethanol_bg1(:,:)
            c_nonsoluble_bg1_pfts(:,:,ipft)  = c_nonsoluble_bg1(:,:)
            c_humus_1_pfts(:,:,ipft)         = c_humus_1(:,:)
            ! size class 2
            c_acid_ag2_pfts(:,:,ipft)        = c_acid_ag2(:,:)
            c_water_ag2_pfts(:,:,ipft)       = c_water_ag2(:,:)
            c_ethanol_ag2_pfts(:,:,ipft)     = c_ethanol_ag2(:,:)
            c_nonsoluble_ag2_pfts(:,:,ipft)  = c_nonsoluble_ag2(:,:)
            c_acid_bg2_pfts(:,:,ipft)        = c_acid_bg2(:,:)
            c_water_bg2_pfts(:,:,ipft)       = c_water_bg2(:,:)
            c_ethanol_bg2_pfts(:,:,ipft)     = c_ethanol_bg2(:,:)
            c_nonsoluble_bg2_pfts(:,:,ipft)  = c_nonsoluble_bg2(:,:)
            c_humus_2_pfts(:,:,ipft)         = c_humus_2(:,:)
            c_bg_sum_pfts(:,:,ipft)          = c_bg_sum(:,:)
            c_ag_sum_1_pfts(:,:,ipft)        = c_ag_sum_1(:,:)
            c_litter_green_ag_pfts(:,:,ipft) = c_litter_green_ag(:,:)
            c_litter_green_bg_pfts(:,:,ipft) = c_litter_green_bg(:,:)
            c_litter_wood_ag_pfts(:,:,ipft)  = c_litter_wood_ag(:,:)
            c_litter_wood_bg_pfts(:,:,ipft)  = c_litter_wood_bg(:,:)
            !
            NPP_pot_sum_pfts(:,:,ipft)          = NPP_pot_sum(:,:)
            GPP_sum_pfts(:,:,ipft)              = GPP_sum(:,:)
            alpha_pfts(:,:,ipft)                = alpha(:,:)
            alpha_sum_pfts(:,:,ipft)            = alpha_sum(:,:)
            alpha_yDayMean_pfts(:,:,ipft)       = alpha_yDayMean(:,:)
            soil_respiration_pfts(:,:,ipft)     = soil_respiration(:,:)
            NPP_flux_correction_pfts(:,:,ipft)  = NPP_flux_correction(:,:)
            NPP_pot_yDayMean_pfts(:,:,ipft)     = NPP_pot_yDayMean(:,:)
            NPP_act_yDayMean_pfts(:,:,ipft)     = NPP_act_yDayMean(:,:)
            excess_NPP_pfts(:,:,ipft)           = excess_NPP(:,:)
            GPP_yDayMean_pfts(:,:,ipft)         = GPP_yDayMean(:,:)
            root_exudates_pfts(:,:,ipft)        = root_exudates(:,:)
            cflux_c_green_2_herb_pfts(:,:,ipft)           = cflux_c_green_2_herb(:,:)
            cflux_herb_2_littergreen_pfts(:,:,ipft)       = cflux_herb_2_littergreen(:,:)
            cflux_herb_2_atm_pfts(:,:,ipft)               = cflux_herb_2_atm(:,:)
            cflux_dist_greenreserve_2_soil_pfts(:,:,ipft) = cflux_dist_greenreserve_2_soil(:,:)
            cflux_dist_woods_2_soil_pfts(:,:,ipft)        = cflux_dist_woods_2_soil(:,:)
            cflux_fire_all_2_atm_pfts(:,:,ipft)           = cflux_fire_all_2_atm(:,:)
            IF ( dsl4jsb_Config(PHENO_)%active == .TRUE.  .AND.   dsl4jsb_Config(PHENO_)%l_forestRegrowth ) THEN
              current_max_green_pfts(:,:,ipft)            = current_max_green(:,:)
              veg_carbon_at_max_green_pfts(:,:,ipft)      = veg_carbon_at_max_green(:,:)
            END IF

          END DO ! ipft

          ! Calculate the average over all lost areas
          dummy(:,:) = SUM(c_green_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)           / fract_lost_sum(:,:)
          c_green_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_woods_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)           / fract_lost_sum(:,:)
          c_woods_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_reserve_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)         / fract_lost_sum(:,:)
          c_reserve_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_crop_harvest_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)    / fract_lost_sum(:,:)
          c_crop_harvest_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_acid_ag1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
          c_acid_ag1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_water_ag1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          c_water_ag1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_ethanol_ag1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          c_ethanol_ag1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_nonsoluble_ag1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          c_nonsoluble_ag1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_acid_bg1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
          c_acid_bg1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_water_bg1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          c_water_bg1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_ethanol_bg1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          c_ethanol_bg1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_nonsoluble_bg1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          c_nonsoluble_bg1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_humus_1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)         / fract_lost_sum(:,:)
          c_humus_1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_acid_ag2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
          c_acid_ag2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_water_ag2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          c_water_ag2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_water_ag2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          c_water_ag2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_ethanol_ag2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          c_ethanol_ag2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_nonsoluble_ag2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          c_nonsoluble_ag2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_acid_bg2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
          c_acid_bg2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_water_bg2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          c_water_bg2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_ethanol_bg2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          c_ethanol_bg2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_nonsoluble_bg2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          c_nonsoluble_bg2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_humus_2_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)         / fract_lost_sum(:,:)
          c_humus_2_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_bg_sum_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)          / fract_lost_sum(:,:)
          c_bg_sum_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_ag_sum_1_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
          c_ag_sum_1_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_litter_green_ag_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          c_litter_green_ag_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_litter_green_bg_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          c_litter_green_bg_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_litter_wood_ag_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          c_litter_wood_ag_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(c_litter_wood_bg_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          c_litter_wood_bg_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !
          dummy(:,:) = SUM(NPP_pot_sum_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          NPP_pot_sum_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(GPP_sum_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)           / fract_lost_sum(:,:)
          GPP_sum_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(alpha_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)             / fract_lost_sum(:,:)
          alpha_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(alpha_sum_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)         / fract_lost_sum(:,:)
          alpha_sum_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(alpha_yDayMean_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)    / fract_lost_sum(:,:)
          alpha_yDayMean_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(soil_respiration_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          soil_respiration_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(NPP_flux_correction_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          NPP_flux_correction_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(NPP_pot_yDayMean_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          NPP_pot_yDayMean_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(NPP_act_yDayMean_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          NPP_act_yDayMean_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(excess_NPP_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
          excess_NPP_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(GPP_yDayMean_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)      / fract_lost_sum(:,:)
          GPP_yDayMean_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(root_exudates_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          root_exudates_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(cflux_c_green_2_herb_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          cflux_c_green_2_herb_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(cflux_herb_2_littergreen_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          cflux_herb_2_littergreen_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(cflux_herb_2_atm_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)         / fract_lost_sum(:,:)
          cflux_herb_2_atm_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(cflux_dist_greenreserve_2_soil_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          cflux_dist_greenreserve_2_soil_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(cflux_dist_woods_2_soil_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          cflux_dist_woods_2_soil_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(cflux_fire_all_2_atm_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          cflux_fire_all_2_atm_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          IF ( dsl4jsb_Config(PHENO_)%active == .TRUE.  .AND.   dsl4jsb_Config(PHENO_)%l_forestRegrowth ) THEN
            dummy(:,:) = SUM(current_max_green_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
            current_max_green_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(veg_carbon_at_max_green_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
            veg_carbon_at_max_green_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          END IF

          ! Only for carbon conservation check
          C_pools_sum_before(:,:) =                                                             &
            & SUM( fract_pfts_before_dynvegstep * (                                             &
            & c_green_pfts    + c_woods_pfts     + c_reserve_pfts     + c_crop_harvest_pfts   + &
            & c_acid_ag1_pfts + c_water_ag1_pfts + c_ethanol_ag1_pfts + c_nonsoluble_ag1_pfts + &
            & c_acid_ag2_pfts + c_water_ag2_pfts + c_ethanol_ag2_pfts + c_nonsoluble_ag2_pfts + &
            & c_acid_bg1_pfts + c_water_bg1_pfts + c_ethanol_bg1_pfts + c_nonsoluble_bg1_pfts + &
            & c_acid_bg2_pfts + c_water_bg2_pfts + c_ethanol_bg2_pfts + c_nonsoluble_bg2_pfts + &
            & c_humus_1_pfts  + c_humus_2_pfts )                                                &
            & , DIM=3)

          ! Now we can change the variables on all pfts
          CALL relocate_vars_to_new_pft_fractions(   &
            & d_fract_grow_pfts(:,:,:),              & ! in
            & f_old_to_new_grow_pfts(:,:,:),         &
            & f_d_to_new_grow_pfts(:,:,:),           &
            ! avgs
            & c_green_lost_avg                        = c_green_lost_avg(:,:,:),             & ! in
            & c_woods_lost_avg                        = c_woods_lost_avg(:,:,:),             &
            & c_reserve_lost_avg                      = c_reserve_lost_avg(:,:,:),           &
            & c_crop_harvest_lost_avg                 = c_crop_harvest_lost_avg(:,:,:),      &
            ! size class 1
            & c_acid_ag1_lost_avg                     = c_acid_ag1_lost_avg(:,:,:),          &
            & c_water_ag1_lost_avg                    = c_water_ag1_lost_avg(:,:,:),         &
            & c_ethanol_ag1_lost_avg                  = c_ethanol_ag1_lost_avg(:,:,:),       &
            & c_nonsoluble_ag1_lost_avg               = c_nonsoluble_ag1_lost_avg(:,:,:),    &
            & c_acid_bg1_lost_avg                     = c_acid_bg1_lost_avg(:,:,:),          &
            & c_water_bg1_lost_avg                    = c_water_bg1_lost_avg(:,:,:),         &
            & c_ethanol_bg1_lost_avg                  = c_ethanol_bg1_lost_avg(:,:,:),       &
            & c_nonsoluble_bg1_lost_avg               = c_nonsoluble_bg1_lost_avg(:,:,:),    &
            & c_humus_1_lost_avg                      = c_humus_1_lost_avg(:,:,:),           &
            ! size class 2
            & c_acid_ag2_lost_avg                     = c_acid_ag2_lost_avg(:,:,:),          &
            & c_water_ag2_lost_avg                    = c_water_ag2_lost_avg(:,:,:),         &
            & c_ethanol_ag2_lost_avg                  = c_ethanol_ag2_lost_avg(:,:,:),       &
            & c_nonsoluble_ag2_lost_avg               = c_nonsoluble_ag2_lost_avg(:,:,:),    &
            & c_acid_bg2_lost_avg                     = c_acid_bg2_lost_avg(:,:,:),          &
            & c_water_bg2_lost_avg                    = c_water_bg2_lost_avg(:,:,:),         &
            & c_ethanol_bg2_lost_avg                  = c_ethanol_bg2_lost_avg(:,:,:),       &
            & c_nonsoluble_bg2_lost_avg               = c_nonsoluble_bg2_lost_avg(:,:,:),    &
            & c_humus_2_lost_avg                      = c_humus_2_lost_avg(:,:,:),           &
            & c_bg_sum_lost_avg                       = c_bg_sum_lost_avg(:,:,:),            &
            & c_ag_sum_1_lost_avg                     = c_ag_sum_1_lost_avg(:,:,:),          &
            & c_litter_green_ag_lost_avg              = c_litter_green_ag_lost_avg(:,:,:),   &
            & c_litter_green_bg_lost_avg              = c_litter_green_bg_lost_avg(:,:,:),   &
            & c_litter_wood_ag_lost_avg               = c_litter_wood_ag_lost_avg(:,:,:),    &
            & c_litter_wood_bg_lost_avg               = c_litter_wood_bg_lost_avg(:,:,:),    &
            !
            & NPP_pot_sum_lost_avg                    = NPP_pot_sum_lost_avg(:,:,:),         &
            & GPP_sum_lost_avg                        = GPP_sum_lost_avg(:,:,:),             &
            & alpha_lost_avg                          = alpha_lost_avg(:,:,:),               &
            & alpha_sum_lost_avg                      = alpha_sum_lost_avg(:,:,:),           &
            & alpha_yDayMean_lost_avg                 = alpha_yDayMean_lost_avg(:,:,:),      &
            & soil_respiration_lost_avg               = soil_respiration_lost_avg(:,:,:),    &
            & NPP_flux_correction_lost_avg            = NPP_flux_correction_lost_avg(:,:,:), &
            & NPP_pot_yDayMean_lost_avg               = NPP_pot_yDayMean_lost_avg(:,:,:),    &
            & NPP_act_yDayMean_lost_avg               = NPP_act_yDayMean_lost_avg(:,:,:),    &
            & excess_NPP_lost_avg                     = excess_NPP_lost_avg(:,:,:),          &
            & GPP_yDayMean_lost_avg                   = GPP_yDayMean_lost_avg(:,:,:),        &
            & root_exudates_lost_avg                  = root_exudates_lost_avg(:,:,:),       &
            & cflux_c_green_2_herb_lost_avg           = cflux_c_green_2_herb_lost_avg(:,:,:),           &
            & cflux_herb_2_littergreen_lost_avg       = cflux_herb_2_littergreen_lost_avg(:,:,:),       &
            & cflux_herb_2_atm_lost_avg               = cflux_herb_2_atm_lost_avg(:,:,:),               &
            & cflux_dist_greenreserve_2_soil_lost_avg = cflux_dist_greenreserve_2_soil_lost_avg(:,:,:), &
            & cflux_dist_woods_2_soil_lost_avg        = cflux_dist_woods_2_soil_lost_avg(:,:,:),        &
            & cflux_fire_all_2_atm_lost_avg           = cflux_fire_all_2_atm_lost_avg(:,:,:),           &
            & current_max_green_lost_avg              = current_max_green_lost_avg(:,:,:),              &
            & veg_carbon_at_max_green_lost_avg        = veg_carbon_at_max_green_lost_avg(:,:,:),        &
            ! pfts
            & c_green_ = c_green_pfts(:,:,:), c_woods_ = c_woods_pfts(:,:,:), c_reserve_ = c_reserve_pfts(:,:,:),     & ! inout
            & c_crop_harvest_ = c_crop_harvest_pfts(:,:,:),                                                           &
            & c_acid_ag1_ = c_acid_ag1_pfts(:,:,:), c_water_ag1_ = c_water_ag1_pfts(:,:,:),                           &
            & c_ethanol_ag1_ = c_ethanol_ag1_pfts(:,:,:), c_nonsoluble_ag1_ = c_nonsoluble_ag1_pfts(:,:,:),           &
            & c_acid_ag2_ = c_acid_ag2_pfts(:,:,:), c_water_ag2_ = c_water_ag2_pfts(:,:,:),                           &
            & c_ethanol_ag2_ = c_ethanol_ag2_pfts(:,:,:), c_nonsoluble_ag2_ = c_nonsoluble_ag2_pfts(:,:,:),           &
            & c_acid_bg1_ = c_acid_bg1_pfts(:,:,:), c_water_bg1_ = c_water_bg1_pfts(:,:,:),                           &
            & c_ethanol_bg1_ = c_ethanol_bg1_pfts(:,:,:), c_nonsoluble_bg1_ = c_nonsoluble_bg1_pfts(:,:,:),           &
            & c_acid_bg2_ = c_acid_bg2_pfts(:,:,:), c_water_bg2_ = c_water_bg2_pfts(:,:,:),                           &
            & c_ethanol_bg2_ = c_ethanol_bg2_pfts(:,:,:), c_nonsoluble_bg2_ = c_nonsoluble_bg2_pfts(:,:,:),           &
            & c_humus_1_ = c_humus_1_pfts(:,:,:), c_humus_2_ = c_humus_2_pfts(:,:,:),                                 &
            & c_bg_sum_ = c_bg_sum_pfts(:,:,:), c_ag_sum_1_ = c_ag_sum_1_pfts(:,:,:),                                 &
            & c_litter_green_ag_ = c_litter_green_ag_pfts(:,:,:), c_litter_green_bg_ = c_litter_green_bg_pfts(:,:,:), &
            & c_litter_wood_ag_ =  c_litter_wood_ag_pfts(:,:,:), c_litter_wood_bg_ = c_litter_wood_bg_pfts(:,:,:),    &
            !
            & NPP_pot_sum_         = NPP_pot_sum_pfts(:,:,:),               GPP_sum_    = GPP_sum_pfts(:,:,:),               &
            & alpha_ = alpha_pfts(:,:,:),  alpha_sum_ = alpha_sum_pfts(:,:,:), alpha_yDayMean_ = alpha_yDayMean_pfts(:,:,:), &
            & soil_respiration_    = soil_respiration_pfts(:,:,:),                                                           &
            & NPP_flux_correction_ = NPP_flux_correction_pfts(:,:,:), NPP_pot_yDayMean_ = NPP_pot_yDayMean_pfts(:,:,:),      &
            & NPP_act_yDayMean_    = NPP_act_yDayMean_pfts(:,:,:),          excess_NPP_ = excess_NPP_pfts(:,:,:),            &
            & GPP_yDayMean_                   = GPP_yDayMean_pfts(:,:,:),                                                    &
            & root_exudates_                  = root_exudates_pfts(:,:,:),                                                   &
            & cflux_c_green_2_herb_           = cflux_c_green_2_herb_pfts(:,:,:),                                            &
            & cflux_herb_2_littergreen_       = cflux_herb_2_littergreen_pfts(:,:,:),                                        &
            & cflux_herb_2_atm_               = cflux_herb_2_atm_pfts(:,:,:),                                                &
            & cflux_dist_greenreserve_2_soil_ = cflux_dist_greenreserve_2_soil_pfts(:,:,:),                                  &
            & cflux_dist_woods_2_soil_        = cflux_dist_woods_2_soil_pfts(:,:,:),                                         &
            & cflux_fire_all_2_atm_           = cflux_fire_all_2_atm_pfts(:,:,:),                                            &
            & current_max_green_              = current_max_green_pfts(:,:,:),                                               &
            & veg_carbon_at_max_green_        = veg_carbon_at_max_green_pfts(:,:,:)                                          &
            & )


          ! Carbon conservation check
          C_pools_sum_after(:,:) =                                                              &
            & SUM( fract_pfts_currentyear_maps   * (                                            &
            & c_green_pfts    + c_woods_pfts     + c_reserve_pfts     + c_crop_harvest_pfts   + &
            & c_acid_ag1_pfts + c_water_ag1_pfts + c_ethanol_ag1_pfts + c_nonsoluble_ag1_pfts + &
            & c_acid_ag2_pfts + c_water_ag2_pfts + c_ethanol_ag2_pfts + c_nonsoluble_ag2_pfts + &
            & c_acid_bg1_pfts + c_water_bg1_pfts + c_ethanol_bg1_pfts + c_nonsoluble_bg1_pfts + &
            & c_acid_bg2_pfts + c_water_bg2_pfts + c_ethanol_bg2_pfts + c_nonsoluble_bg2_pfts + &
            & c_humus_1_pfts  + c_humus_2_pfts )                                                &
            & , DIM=3)

          dynveg__mem%C_pools_sum_check%ptr(:,:) = C_pools_sum_after(:,:) - C_pools_sum_before(:,:)
          ! DO i=1,hgrid%nproma
          !   DO j=1,hgrid%nblks
          !       !IF ( ABS(C_pools_sum_check(i,j))  > EPSILON(1._wp) ) THEN
          !       !  write(*,*) 'C_pools_sum_check of dynveg maps > or < EPSILON(1._wp):' , C_pools_sum_check(i,j)
          !       !END IF
          !       IF ( ABS(dynveg__mem%C_pools_sum_check%ptr(i,j))  > 2*EPSILON(1._wp) ) THEN
          !         write(*,*) 'C_pools_sum_check of dynveg maps > or < EPSILON(1._wp):' , dynveg__mem%C_pools_sum_check%ptr(i,j)
          !       END IF
          !   END DO ! j
          ! END DO ! i
          ! IF (debug_on()) CALL message(TRIM(routine), '...C_pools_sum_check of dynveg maps corresponding to '//TRIM(filename_maps)//'finished...')


          ! Hand down the results to child tiles
          ! Also possible with mo_jsb_memory_class.f90: SUBROUTINE HandDown_var_real3d
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(CARBON_)%p)
                  CLASS IS (t_carbon_memory)
                    carbon__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_carbon_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            IF (iarea == 1 ) THEN ! _ca variables
              carbon__mem%c_green%ptr(:,:)   = c_green_pfts(:,:,ipft)
              carbon__mem%c_reserve%ptr(:,:) = c_reserve_pfts(:,:,ipft)
              carbon__mem%c_woods%ptr(:,:)   = c_woods_pfts(:,:,ipft)
              carbon__mem%c_crop_harvest%ptr(:,:)    = c_crop_harvest_pfts(:,:,ipft)
              carbon__mem%c_litter_green_ag%ptr(:,:) = c_litter_green_ag_pfts(:,:,ipft)
              carbon__mem%c_litter_green_bg%ptr(:,:) = c_litter_green_bg_pfts(:,:,ipft)
              carbon__mem%c_litter_wood_ag %ptr(:,:) = c_litter_wood_ag_pfts(:,:,ipft)
              carbon__mem%c_litter_wood_bg %ptr(:,:) = c_litter_wood_bg_pfts(:,:,ipft)
              ! size class 1
              carbon__mem%c_acid_ag1%ptr(:,:)        = c_acid_ag1_pfts(:,:,ipft)
              carbon__mem%c_water_ag1%ptr(:,:)       = c_water_ag1_pfts(:,:,ipft)
              carbon__mem%c_ethanol_ag1%ptr(:,:)     = c_ethanol_ag1_pfts(:,:,ipft)
              carbon__mem%c_nonsoluble_ag1%ptr(:,:)  = c_nonsoluble_ag1_pfts(:,:,ipft)
              carbon__mem%c_acid_bg1%ptr(:,:)        = c_acid_bg1_pfts(:,:,ipft)
              carbon__mem%c_water_bg1%ptr(:,:)       = c_water_bg1_pfts(:,:,ipft)
              carbon__mem%c_ethanol_bg1%ptr(:,:)     = c_ethanol_bg1_pfts(:,:,ipft)
              carbon__mem%c_nonsoluble_bg1%ptr(:,:)  = c_nonsoluble_bg1_pfts(:,:,ipft)
              carbon__mem%c_humus_1%ptr(:,:)         = c_humus_1_pfts(:,:,ipft)
              ! size class 2
              carbon__mem%c_acid_ag2%ptr(:,:)        = c_acid_ag2_pfts(:,:,ipft)
              carbon__mem%c_water_ag2%ptr(:,:)       = c_water_ag2_pfts(:,:,ipft)
              carbon__mem%c_ethanol_ag2%ptr(:,:)     = c_ethanol_ag2_pfts(:,:,ipft)
              carbon__mem%c_nonsoluble_ag2%ptr(:,:)  = c_nonsoluble_ag2_pfts(:,:,ipft)
              carbon__mem%c_acid_bg2%ptr(:,:)        = c_acid_bg2_pfts(:,:,ipft)
              carbon__mem%c_water_bg2%ptr(:,:)       = c_water_bg2_pfts(:,:,ipft)
              carbon__mem%c_ethanol_bg2%ptr(:,:)     = c_ethanol_bg2_pfts(:,:,ipft)
              carbon__mem%c_nonsoluble_bg2%ptr(:,:)  = c_nonsoluble_bg2_pfts(:,:,ipft)
              carbon__mem%c_humus_2%ptr(:,:)         = c_humus_2_pfts(:,:,ipft)
              carbon__mem%c_bg_sum%ptr(:,:)          = c_bg_sum_pfts(:,:,ipft)
              carbon__mem%c_ag_sum_1%ptr(:,:)        = c_ag_sum_1_pfts(:,:,ipft)
              carbon__mem%c_litter_green_ag%ptr(:,:) = c_litter_green_ag_pfts(:,:,ipft)
              carbon__mem%c_litter_green_bg%ptr(:,:) = c_litter_green_bg_pfts(:,:,ipft)
              carbon__mem%c_litter_wood_ag%ptr(:,:)  = c_litter_wood_ag_pfts(:,:,ipft)
              carbon__mem%c_litter_wood_bg%ptr(:,:)  = c_litter_wood_bg_pfts(:,:,ipft)
              !
              carbon__mem%NPP_pot_sum%ptr(:,:)          = NPP_pot_sum_pfts(:,:,ipft)
              carbon__mem%GPP_sum%ptr(:,:)              = GPP_sum_pfts(:,:,ipft)
              carbon__mem%alpha%ptr(:,:)                = alpha_pfts(:,:,ipft)
              carbon__mem%alpha_sum%ptr(:,:)            = alpha_sum_pfts(:,:,ipft)
              carbon__mem%alpha_yDayMean%ptr(:,:)       = alpha_yDayMean_pfts(:,:,ipft)
              carbon__mem%soil_respiration%ptr(:,:)     = soil_respiration_pfts(:,:,ipft)
              carbon__mem%NPP_flux_correction%ptr(:,:)  = NPP_flux_correction_pfts(:,:,ipft)
              carbon__mem%NPP_pot_yDayMean%ptr(:,:)     = NPP_pot_yDayMean_pfts(:,:,ipft)
              carbon__mem%NPP_act_yDayMean%ptr(:,:)     = NPP_act_yDayMean_pfts(:,:,ipft)
              carbon__mem%excess_NPP%ptr(:,:)           = excess_NPP_pfts(:,:,ipft)
              carbon__mem%GPP_yDayMean%ptr(:,:)         = GPP_yDayMean_pfts(:,:,ipft)
              carbon__mem%root_exudates%ptr(:,:)        = root_exudates_pfts(:,:,ipft)
              carbon__mem%cflux_c_green_2_herb%ptr(:,:) = cflux_c_green_2_herb_pfts(:,:,ipft)
              carbon__mem%cflux_herb_2_littergreen%ptr(:,:) = cflux_herb_2_littergreen_pfts(:,:,ipft)
              carbon__mem%cflux_herb_2_atm%ptr(:,:)         = cflux_herb_2_atm_pfts(:,:,ipft)
              carbon__mem%cflux_dist_greenreserve_2_soil%ptr(:,:) = cflux_dist_greenreserve_2_soil_pfts(:,:,ipft)
              carbon__mem%cflux_dist_woods_2_soil%ptr(:,:)        = cflux_dist_woods_2_soil_pfts(:,:,ipft)
              carbon__mem%cflux_fire_all_2_atm%ptr(:,:)           = cflux_fire_all_2_atm_pfts(:,:,ipft)
              IF ( dsl4jsb_Config(PHENO_)%active == .TRUE.  .AND.   dsl4jsb_Config(PHENO_)%l_forestRegrowth ) THEN
                carbon__mem%current_max_green%ptr(:,:)       = current_max_green_pfts(:,:,ipft)
                carbon__mem%veg_carbon_at_max_green%ptr(:,:) = veg_carbon_at_max_green_pfts(:,:,ipft)
              END IF
            ELSE ! _ta variables
              carbon__mem%c_green_ta%ptr(:,:)           = c_green_pfts(:,:,ipft)
              carbon__mem%c_reserve_ta%ptr(:,:)         = c_reserve_pfts(:,:,ipft)
              carbon__mem%c_woods_ta%ptr(:,:)           = c_woods_pfts(:,:,ipft)
              carbon__mem%c_crop_harvest_ta%ptr(:,:)    = c_crop_harvest_pfts(:,:,ipft)
              carbon__mem%c_litter_green_ag_ta%ptr(:,:) = c_litter_green_ag_pfts(:,:,ipft)
              carbon__mem%c_litter_green_bg_ta%ptr(:,:) = c_litter_green_bg_pfts(:,:,ipft)
              carbon__mem%c_litter_wood_ag_ta%ptr(:,:)  = c_litter_wood_ag_pfts(:,:,ipft)
              carbon__mem%c_litter_wood_bg_ta%ptr(:,:)  = c_litter_wood_bg_pfts(:,:,ipft)
              !SIZE CLASS 1
              carbon__mem%c_acid_ag1_ta%ptr(:,:)        = c_acid_ag1_pfts(:,:,ipft)
              carbon__mem%c_water_ag1_ta%ptr(:,:)       = c_water_ag1_pfts(:,:,ipft)
              carbon__mem%c_ethanol_ag1_ta%ptr(:,:)     = c_ethanol_ag1_pfts(:,:,ipft)
              carbon__mem%c_nonsoluble_ag1_ta%ptr(:,:)  = c_nonsoluble_ag1_pfts(:,:,ipft)
              carbon__mem%c_acid_bg1_ta%ptr(:,:)        = c_acid_bg1_pfts(:,:,ipft)
              carbon__mem%c_water_bg1_ta%ptr(:,:)       = c_water_bg1_pfts(:,:,ipft)
              carbon__mem%c_ethanol_bg1_ta%ptr(:,:)     = c_ethanol_bg1_pfts(:,:,ipft)
              carbon__mem%c_nonsoluble_bg1_ta%ptr(:,:)  = c_nonsoluble_bg1_pfts(:,:,ipft)
              carbon__mem%c_humus_1_ta%ptr(:,:)         = c_humus_1_pfts(:,:,ipft)
              !SIZE CLASS 2
              carbon__mem%c_acid_ag2_ta%ptr(:,:)        = c_acid_ag2_pfts(:,:,ipft)
              carbon__mem%c_water_ag2_ta%ptr(:,:)       = c_water_ag2_pfts(:,:,ipft)
              carbon__mem%c_ethanol_ag2_ta%ptr(:,:)     = c_ethanol_ag2_pfts(:,:,ipft)
              carbon__mem%c_nonsoluble_ag2_ta%ptr(:,:)  = c_nonsoluble_ag2_pfts(:,:,ipft)
              carbon__mem%c_acid_bg2_ta%ptr(:,:)        = c_acid_bg2_pfts(:,:,ipft)
              carbon__mem%c_water_bg2_ta%ptr(:,:)       = c_water_bg2_pfts(:,:,ipft)
              carbon__mem%c_ethanol_bg2_ta%ptr(:,:)     = c_ethanol_bg2_pfts(:,:,ipft)
              carbon__mem%c_nonsoluble_bg2_ta%ptr(:,:)  = c_nonsoluble_bg2_pfts(:,:,ipft)
              carbon__mem%c_humus_2_ta%ptr(:,:)         = c_humus_2_pfts(:,:,ipft)
              ! no _ta var: carbon__mem%c_bg_sum_ta%ptr(:,:)         = c_bg_sum_pfts(:,:,ipft)
              ! no _ta var: carbon__mem%c_ag_sum_1_ta%ptr(:,:)       = c_ag_sum_1_pfts(:,:,ipft)
              carbon__mem%c_litter_green_ag_ta%ptr(:,:) = c_litter_green_ag_pfts(:,:,ipft)
              carbon__mem%c_litter_green_bg_ta%ptr(:,:) = c_litter_green_bg_pfts(:,:,ipft)
              carbon__mem%c_litter_wood_ag_ta%ptr(:,:)  = c_litter_wood_ag_pfts(:,:,ipft)
              carbon__mem%c_litter_wood_bg_ta%ptr(:,:)  = c_litter_wood_bg_pfts(:,:,ipft)
              !
              ! no _ta var: carbon__mem%NPP_pot_sum_ta%ptr(:,:)          = NPP_pot_sum_pfts(:,:,ipft)
              ! no _ta var: carbon__mem%GPP_sum_ta%ptr(:,:)              = GPP_sum_pfts(:,:,ipft)
              ! no _ta var: carbon__mem%alpha_ta%ptr(:,:)                = alpha_pfts(:,:,ipft)
              ! no _ta var: carbon__mem%alpha_sum_ta%ptr(:,:)            = alpha_sum_pfts(:,:,ipft)
              ! no _ta var: carbon__mem%alpha_yDayMean_ta%ptr(:,:)       = alpha_yDayMean_pfts(:,:,ipft)
              carbon__mem%soil_respiration_ta%ptr(:,:)     = soil_respiration_pfts(:,:,ipft)
              carbon__mem%NPP_flux_correction_ta%ptr(:,:)  = NPP_flux_correction_pfts(:,:,ipft)
              carbon__mem%NPP_pot_yDayMean_ta%ptr(:,:)     = NPP_pot_yDayMean_pfts(:,:,ipft)
              carbon__mem%NPP_act_yDayMean_ta%ptr(:,:)     = NPP_act_yDayMean_pfts(:,:,ipft)
              carbon__mem%excess_NPP_ta%ptr(:,:)           = excess_NPP_pfts(:,:,ipft)
              carbon__mem%GPP_yDayMean_ta%ptr(:,:)         = GPP_yDayMean_pfts(:,:,ipft)
              carbon__mem%root_exudates_ta%ptr(:,:)        = root_exudates_pfts(:,:,ipft)
              carbon__mem%cflux_c_green_2_herb_ta%ptr(:,:) = cflux_c_green_2_herb_pfts(:,:,ipft)
              ! no _ta var: carbon__mem%cflux_herb_2_littergreen_ta%ptr(:,:) = cflux_herb_2_littergreen_pfts(:,:,ipft)
              ! no _ta var: carbon__mem%cflux_herb_2_atm_ta%ptr(:,:)         = cflux_herb_2_atm_pfts(:,:,ipft)
              carbon__mem%cflux_dist_greenreserve_2_soil_ta%ptr(:,:) = cflux_dist_greenreserve_2_soil_pfts(:,:,ipft)
              carbon__mem%cflux_dist_woods_2_soil_ta%ptr(:,:)        = cflux_dist_woods_2_soil_pfts(:,:,ipft)
              carbon__mem%cflux_fire_all_2_atm_ta%ptr(:,:)           = cflux_fire_all_2_atm_pfts(:,:,ipft)
              !IF ( dsl4jsb_Config(PHENO_)%active == .TRUE.  .AND.   dsl4jsb_Config(PHENO_)%l_forestRegrowth ) THEN
                ! no _ta var: carbon__mem%current_max_green_ta%ptr(:,:)       = current_max_green_pfts(:,:,ipft)
                ! no _ta var: carbon__mem%veg_carbon_at_max_green_ta%ptr(:,:) = veg_carbon_at_max_green_pfts(:,:,ipft)
              !END IF
            END IF

          END DO ! ipft
        END DO ! iarea
      END IF ! CARBON is active


      IF (dsl4jsb_Config(ASSIMI_)%active == .TRUE.) THEN
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(ASSIMI_)%p)
                  CLASS IS (t_assimi_memory)
                    assimi__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_assimi_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            gross_assimilation_ca   => assimi__mem%gross_assimilation_ca%ptr(:,:) ! inout
            gross_assimilation      => assimi__mem%gross_assimilation%ptr(:,:)    ! inout
            dark_respiration_ca     => assimi__mem%dark_respiration_ca%ptr(:,:)   ! inout
            dark_respiration        => assimi__mem%dark_respiration%ptr(:,:)      ! inout
            NPP_pot_rate_ca         => assimi__mem%NPP_pot_rate_ca%ptr(:,:)       ! inout
            canopy_cond             => assimi__mem%canopy_cond_unlimited%ptr(:,:) ! inout
            cover_fract_pot         => assimi__mem%cover_fract_pot%ptr(:,:)       ! inout

            assimi_scaling_fact_cl  => assimi__mem%assimi_scaling_fact_cl%ptr(:,:,:)  ! inout
            canopy_cond_cl          => assimi__mem%canopy_cond_cl%ptr(:,:,:)          ! inout
            gross_assimilation_cl   => assimi__mem%gross_assimilation_cl%ptr(:,:,:)   ! inout
            dark_respiration_cl     => assimi__mem%dark_respiration_cl%ptr(:,:,:)     ! inout
            carbox_rate_max_cl      => assimi__mem%carbox_rate_max_cl%ptr(:,:,:)      ! inout
            e_transport_rate_max_cl => assimi__mem%e_transport_rate_max_cl%ptr(:,:,:) ! inout
            carbox_rate_cl          => assimi__mem%carbox_rate_cl%ptr(:,:,:)          ! inout
            e_transport_rate_cl     => assimi__mem%e_transport_rate_cl%ptr(:,:,:)     ! inout
            CO2_conc_leaf_cl        => assimi__mem%CO2_conc_leaf_cl%ptr(:,:,:)        ! inout

            ! Collect
            gross_assimilation_ca_pfts(:,:,ipft)     = gross_assimilation_ca(:,:)
            gross_assimilation_pfts(:,:,ipft)        = gross_assimilation(:,:)
            dark_respiration_ca_pfts(:,:,ipft)       = dark_respiration_ca(:,:)
            dark_respiration_pfts(:,:,ipft)          = dark_respiration(:,:)
            NPP_pot_rate_ca_pfts(:,:,ipft)           = NPP_pot_rate_ca(:,:)
            canopy_cond_pfts(:,:,ipft)     = canopy_cond(:,:)
            cover_fract_pot_pfts(:,:,ipft)           = cover_fract_pot(:,:)

            assimi_scaling_fact_cl_pfts(:,:,:,ipft)  = assimi_scaling_fact_cl(:,:,:)
            canopy_cond_cl_pfts(:,:,:,ipft)          = canopy_cond_cl(:,:,:)
            gross_assimilation_cl_pfts(:,:,:,ipft)   = gross_assimilation_cl(:,:,:)
            dark_respiration_cl_pfts(:,:,:,ipft)     = dark_respiration_cl(:,:,:)
            carbox_rate_max_cl_pfts(:,:,:,ipft)      = carbox_rate_max_cl(:,:,:)
            e_transport_rate_max_cl_pfts(:,:,:,ipft) = e_transport_rate_max_cl(:,:,:)
            carbox_rate_cl_pfts(:,:,:,ipft)          = carbox_rate_cl(:,:,:)
            e_transport_rate_cl_pfts(:,:,:,ipft)     = e_transport_rate_cl(:,:,:)
            CO2_conc_leaf_cl_pfts(:,:,:,ipft)        = CO2_conc_leaf_cl(:,:,:)
          END DO ! ipft

          ! Calculate the average over all lost areas
          dummy(:,:) = SUM(gross_assimilation_ca_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          gross_assimilation_ca_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(gross_assimilation_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)    / fract_lost_sum(:,:)
          gross_assimilation_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(dark_respiration_ca_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)   / fract_lost_sum(:,:)
          dark_respiration_ca_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(dark_respiration_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)      / fract_lost_sum(:,:)
          dark_respiration_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(NPP_pot_rate_ca_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          NPP_pot_rate_ca_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(canopy_cond_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          canopy_cond_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(cover_fract_pot_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          cover_fract_pot_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )

          DO icanopy=1,ncanopy
            dummy(:,:) = SUM(assimi_scaling_fact_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
            assimi_scaling_fact_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(canopy_cond_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)          / fract_lost_sum(:,:)
            canopy_cond_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(gross_assimilation_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)   / fract_lost_sum(:,:)
            gross_assimilation_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(dark_respiration_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
            dark_respiration_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(carbox_rate_max_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)      / fract_lost_sum(:,:)
            carbox_rate_max_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(e_transport_rate_max_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
            e_transport_rate_max_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(carbox_rate_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)          / fract_lost_sum(:,:)
            carbox_rate_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(e_transport_rate_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
            e_transport_rate_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(CO2_conc_leaf_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
            CO2_conc_leaf_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          END DO

          ! Now we can change the variables on all pfts
          CALL relocate_vars_to_new_pft_fractions( &
            & d_fract_grow_pfts(:,:,:),            &
            & f_old_to_new_grow_pfts(:,:,:),       &
            & f_d_to_new_grow_pfts(:,:,:),         &
            ! avgs
            & gross_assimilation_ca_lost_avg = gross_assimilation_ca_lost_avg(:,:,:), &
            & gross_assimilation_lost_avg = gross_assimilation_lost_avg(:,:,:),       &
            & dark_respiration_ca_lost_avg = dark_respiration_ca_lost_avg(:,:,:),     &
            & dark_respiration_lost_avg = dark_respiration_lost_avg(:,:,:),           &
            & NPP_pot_rate_ca_lost_avg = NPP_pot_rate_ca_lost_avg(:,:,:),             &
            & canopy_cond_lost_avg = canopy_cond_lost_avg(:,:,:),                     &
            & cover_fract_pot_lost_avg = cover_fract_pot_lost_avg(:,:,:),             &
            ! pfts
            & gross_assimilation_ca = gross_assimilation_ca_pfts(:,:,:), &
            & gross_assimilation = gross_assimilation_pfts(:,:,:),       &
            & dark_respiration_ca = dark_respiration_ca_pfts(:,:,:),     &
            & dark_respiration = dark_respiration_pfts(:,:,:),           &
            & NPP_pot_rate_ca = NPP_pot_rate_ca_pfts(:,:,:),             &
            & canopy_cond = canopy_cond_pfts(:,:,:),                     &
            & cover_fract_pot = cover_fract_pot_pfts(:,:,:)              &
            & )

          DO icanopy=1,ncanopy
            CALL relocate_vars_to_new_pft_fractions( &
              & d_fract_grow_pfts(:,:,:),            &
              & f_old_to_new_grow_pfts(:,:,:),       &
              & f_d_to_new_grow_pfts(:,:,:),         &
              ! avgs
              & assimi_scaling_fact_cl_lost_avg = assimi_scaling_fact_cl_lost_avg(:,icanopy,:,:),   &
              & canopy_cond_cl_lost_avg = canopy_cond_cl_lost_avg(:,icanopy,:,:),                   &
              & gross_assimilation_cl_lost_avg = gross_assimilation_cl_lost_avg(:,icanopy,:,:),     &
              & dark_respiration_cl_lost_avg = dark_respiration_cl_lost_avg(:,icanopy,:,:),         &
              & carbox_rate_max_cl_lost_avg = carbox_rate_max_cl_lost_avg(:,icanopy,:,:),           &
              & e_transport_rate_max_cl_lost_avg = e_transport_rate_max_cl_lost_avg(:,icanopy,:,:), &
              & carbox_rate_cl_lost_avg = carbox_rate_cl_lost_avg(:,icanopy,:,:),                   &
              & e_transport_rate_cl_lost_avg = e_transport_rate_cl_lost_avg(:,icanopy,:,:),         &
              & CO2_conc_leaf_cl_lost_avg = CO2_conc_leaf_cl_lost_avg(:,icanopy,:,:),               &
              ! pfts
              & assimi_scaling_fact_cl = assimi_scaling_fact_cl_pfts(:,icanopy,:,:),   &
              & canopy_cond_cl = canopy_cond_cl_pfts(:,icanopy,:,:),                   &
              & gross_assimilation_cl = gross_assimilation_cl_pfts(:,icanopy,:,:),     &
              & dark_respiration_cl = dark_respiration_cl_pfts(:,icanopy,:,:),         &
              & carbox_rate_max_cl = carbox_rate_max_cl_pfts(:,icanopy,:,:),           &
              & e_transport_rate_max_cl = e_transport_rate_max_cl_pfts(:,icanopy,:,:), &
              & carbox_rate_cl = carbox_rate_cl_pfts(:,icanopy,:,:),                   &
              & e_transport_rate_cl = e_transport_rate_cl_pfts(:,icanopy,:,:),         &
              & CO2_conc_leaf_cl = CO2_conc_leaf_cl_pfts(:,icanopy,:,:)                &
              & )
          END DO


          ! Hand down the results to child tiles
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(ASSIMI_)%p)
                  CLASS IS (t_assimi_memory)
                    assimi__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_assimi_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            assimi__mem%gross_assimilation_ca%ptr(:,:)     = gross_assimilation_ca_pfts(:,:,ipft)
            assimi__mem%gross_assimilation%ptr(:,:)        = gross_assimilation_pfts(:,:,ipft)
            assimi__mem%dark_respiration_ca%ptr(:,:)       = dark_respiration_ca_pfts(:,:,ipft)
            assimi__mem%dark_respiration%ptr(:,:)          = dark_respiration_pfts(:,:,ipft)
            assimi__mem%NPP_pot_rate_ca%ptr(:,:)           = NPP_pot_rate_ca_pfts(:,:,ipft)
            assimi__mem%canopy_cond_unlimited%ptr(:,:)     = canopy_cond_pfts(:,:,ipft)
            assimi__mem%cover_fract_pot%ptr(:,:)           = cover_fract_pot_pfts(:,:,ipft)

            assimi__mem%assimi_scaling_fact_cl%ptr(:,:,:)  = assimi_scaling_fact_cl_pfts(:,:,:,ipft)
            assimi__mem%canopy_cond_cl%ptr(:,:,:)          = canopy_cond_cl_pfts(:,:,:,ipft)
            assimi__mem%gross_assimilation_cl%ptr(:,:,:)   = gross_assimilation_cl_pfts(:,:,:,ipft)
            assimi__mem%dark_respiration_cl%ptr(:,:,:)     = dark_respiration_cl_pfts(:,:,:,ipft)
            assimi__mem%carbox_rate_max_cl%ptr(:,:,:)      = carbox_rate_max_cl_pfts(:,:,:,ipft)
            assimi__mem%e_transport_rate_max_cl%ptr(:,:,:) = e_transport_rate_max_cl_pfts(:,:,:,ipft)
            assimi__mem%carbox_rate_cl%ptr(:,:,:)          = carbox_rate_cl_pfts(:,:,:,ipft)
            assimi__mem%e_transport_rate_cl%ptr(:,:,:)     = e_transport_rate_cl_pfts(:,:,:,ipft)
            assimi__mem%CO2_conc_leaf_cl%ptr(:,:,:)        = CO2_conc_leaf_cl_pfts(:,:,:,ipft)

          END DO ! ipft
      END IF ! ASSIMI active


      IF (dsl4jsb_Config(DISTURB_)%active == .TRUE.) THEN
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(DISTURB_)%p)
                  CLASS IS (t_disturb_memory)
                    disturb__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_disturb_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            cconservation_fire => disturb__mem%cconservation_fire%ptr(:,:) ! inout
            cconservation_wind => disturb__mem%cconservation_wind%ptr(:,:) ! inout
            burned_fract       => disturb__mem%burned_fract%ptr(:,:)       ! inout
            damaged_fract      => disturb__mem%damaged_fract%ptr(:,:)      ! inout

            ! Collect
            cconservation_fire_pfts(:,:,ipft) = cconservation_fire(:,:)
            cconservation_wind_pfts(:,:,ipft) = cconservation_wind(:,:)
            burned_fract_pfts(:,:,ipft)       = burned_fract(:,:)
            damaged_fract_pfts(:,:,ipft)      = damaged_fract(:,:)
          END DO ! ipft

          ! Calculate the average over all lost areas
          dummy(:,:) = SUM(cconservation_fire_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          cconservation_fire_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(cconservation_wind_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          cconservation_wind_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(burned_fract_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          burned_fract_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(damaged_fract_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)      / fract_lost_sum(:,:)
          damaged_fract_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )

          ! Now we can change the variables on all pfts
          CALL relocate_vars_to_new_pft_fractions( &
            & d_fract_grow_pfts(:,:,:),            &
            & f_old_to_new_grow_pfts(:,:,:),       &
            & f_d_to_new_grow_pfts(:,:,:),         &
            ! avgs
            & cconservation_fire_lost_avg = cconservation_fire_lost_avg(:,:,:), &
            & cconservation_wind_lost_avg = cconservation_wind_lost_avg(:,:,:), &
            & burned_fract_lost_avg       = burned_fract_lost_avg(:,:,:),       &
            & damaged_fract_lost_avg      = damaged_fract_lost_avg(:,:,:),      &
            ! pfts
            & cconservation_fire = cconservation_fire_pfts(:,:,:), &
            & cconservation_wind = cconservation_wind_pfts(:,:,:), &
            & burned_fract       = burned_fract_pfts(:,:,:),       &
            & damaged_fract      = damaged_fract_pfts(:,:,:)       &
            & )


          ! Hand down the results to child tiles
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(DISTURB_)%p)
                  CLASS IS (t_disturb_memory)
                    disturb__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_disturb_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            disturb__mem%cconservation_fire%ptr(:,:) = cconservation_fire_pfts(:,:,ipft)
            disturb__mem%cconservation_wind%ptr(:,:) = cconservation_wind_pfts(:,:,ipft)
            disturb__mem%burned_fract%ptr(:,:)       = burned_fract_pfts(:,:,ipft)
            disturb__mem%damaged_fract%ptr(:,:)      = damaged_fract_pfts(:,:,ipft)
          END DO ! ipft
      END IF ! DISTURB active


      IF (dsl4jsb_Config(PHENO_)%active == .TRUE.) THEN
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(PHENO_)%p)
                  CLASS IS (t_pheno_memory)
                    pheno__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_pheno_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            fract_fpc                    => pheno__mem%fract_fpc%ptr(:,:)                    ! inout
            !fract_fpc_max                => pheno__mem%fract_fpc_max%ptr(:,:)               ! inout
            previous_day_NPP_pot_rate_ca => pheno__mem%previous_day_NPP_pot_rate_ca%ptr(:,:) ! inout
           !IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
           ! cconservation_allom          => pheno__mem%cconservation_allom%ptr(:,:)          ! inout
            day_NPP_sum                  => pheno__mem%day_NPP_sum%ptr(:,:)                  ! inout

            ! Collect
            fract_fpc_pfts(:,:,ipft)                    = fract_fpc(:,:)
            !fract_fpc_max_pfts(:,:,ipft)                = fract_fpc_max(:,:)
            previous_day_NPP_pot_rate_ca_pfts(:,:,ipft) = previous_day_NPP_pot_rate_ca(:,:)
            !IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
            ! cconservation_allom_pfts(:,:,ipft)          = cconservation_allom(:,:)
            day_NPP_sum_pfts(:,:,ipft)                  = day_NPP_sum(:,:)

          END DO ! ipft

          ! Calculate the average over all lost areas
          dummy(:,:) = SUM(fract_fpc_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)                    / fract_lost_sum(:,:)
          fract_fpc_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(fract_fpc_max_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)                / fract_lost_sum(:,:)
          !fract_fpc_max_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(previous_day_NPP_pot_rate_ca_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          previous_day_NPP_pot_rate_ca_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !IF (dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
          ! dummy(:,:) = SUM(cconservation_allom_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)          / fract_lost_sum(:,:)
          ! cconservation_allom_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(day_NPP_sum_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)                  / fract_lost_sum(:,:)
          day_NPP_sum_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )


          ! Now we can change the variables on all pfts
          CALL relocate_vars_to_new_pft_fractions( &
            & d_fract_grow_pfts(:,:,:),            &
            & f_old_to_new_grow_pfts(:,:,:),       &
            & f_d_to_new_grow_pfts(:,:,:),         &
            ! avgs
            & fract_fpc_lost_avg                    = fract_fpc_lost_avg(:,:,:),                    &
            !& fract_fpc_max_lost_avg               = fract_fpc_max_lost_avg(:,:,:),                 &
            & previous_day_NPP_pot_rate_ca_lost_avg = previous_day_NPP_pot_rate_ca_lost_avg(:,:,:), &
            ! & cconservation_allom_lost_avg          = cconservation_allom_lost_avg(:,:,:),          &
            & day_NPP_sum_lost_avg                  = day_NPP_sum_lost_avg(:,:,:),                  &
            ! pfts
            & fract_fpc                    = fract_fpc_pfts(:,:,:),                    &
            !& fract_fpc_max               = fract_fpc_max_pfts(:,:,:),                 &
            & previous_day_NPP_pot_rate_ca = previous_day_NPP_pot_rate_ca_pfts(:,:,:), &
            ! & cconservation_allom          = cconservation_allom_pfts(:,:,:),          &
            & day_NPP_sum                  = day_NPP_sum_pfts(:,:,:)                   &
            & )

          ! Hand down the results to child tiles
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(PHENO_)%p)
                  CLASS IS (t_pheno_memory)
                    pheno__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_pheno_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            pheno__mem%fract_fpc%ptr(:,:)                    = fract_fpc_pfts(:,:,ipft)
            !pheno__mem%fract_fpc_max%ptr(:,:)               = fract_fpc_max_pfts(:,:,ipft)
            pheno__mem%previous_day_NPP_pot_rate_ca%ptr(:,:) = previous_day_NPP_pot_rate_ca_pfts(:,:,ipft)
            ! pheno__mem%cconservation_allom%ptr(:,:)          = cconservation_allom_pfts(:,:,ipft)
            pheno__mem%day_NPP_sum%ptr(:,:)                  = day_NPP_sum_pfts(:,:,ipft)
          END DO ! ipft
      END IF ! PHENO active


      IF (dsl4jsb_Config(RAD_)%active == .TRUE.) THEN
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(RAD_)%p)
                  CLASS IS (t_rad_memory)
                    rad__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_rad_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            !sw_srf_net    => rad__mem%sw_srf_net%ptr(:,:)    ! inout
            !swvis_srf_net => rad__mem%swvis_srf_net%ptr(:,:) ! inout
            !swnir_srf_net => rad__mem%swnir_srf_net%ptr(:,:) ! inout
            !lw_srf_net    => rad__mem%lw_srf_net%ptr(:,:)    ! inout
            !rad_srf_net   => rad__mem%rad_srf_net%ptr(:,:)   ! inout
            !alb_vis       => rad__mem%alb_vis%ptr(:,:)       ! inout
            !alb_nir       => rad__mem%alb_nir%ptr(:,:)       ! inout
            !alb_vis_soil  => rad__mem%alb_vis_soil%ptr(:,:)  ! inout
            !alb_nir_soil  => rad__mem%alb_nir_soil%ptr(:,:)  ! inout
            !alb_vis_lnd   => rad__mem%alb_vis_lnd%ptr(:,:)   ! inout
            !alb_nir_lnd   => rad__mem%alb_nir_lnd%ptr(:,:)   ! inout
            !alb_vis_can   => rad__mem%alb_vis_can%ptr(:,:)   ! inout
            !alb_nir_can   => rad__mem%alb_nir_can%ptr(:,:)   ! inout
            faPAR                 => rad__mem%faPAR%ptr(:,:)         ! inout
            soil_reflectivity_par => rad__mem%soil_reflectivity_par%ptr(:,:) ! inout
            ! apar_nacc             => rad__mem%apar_nacc%ptr(:,:)             ! inout

            faPAR_cl              => rad__mem%faPAR_cl%ptr(:,:,:)             ! inout
            apar_per_lai_cl       => rad__mem%apar_per_lai_cl%ptr(:,:,:)      ! inout
            ! apar_per_lai_cl_nacc  => rad__mem%apar_per_lai_cl_nacc%ptr(:,:,:) ! inout
            ! apar_per_lai_cl_tavg  => rad__mem%apar_per_lai_cl_tavg%ptr(:,:,:) ! inout

            ! Collect
            !sw_srf_net_pfts(:,:,ipft)             = sw_srf_net(:,:)
            !swvis_srf_net_pfts(:,:,ipft)          = swvis_srf_net(:,:)
            !swnir_srf_net_pfts(:,:,ipft)          = swnir_srf_net(:,:)
            !lw_srf_net_pfts(:,:,ipft)             = lw_srf_net(:,:)
            !rad_srf_net_pfts(:,:,ipft)            = rad_srf_net(:,:)
            !alb_vis_pfts(:,:,ipft)                = alb_vis(:,:)
            !alb_nir_pfts(:,:,ipft)                = alb_nir(:,:)
            !alb_vis_soil_pfts(:,:,ipft)           = alb_vis_soil(:,:)
            !alb_nir_soil_pfts(:,:,ipft)           = alb_nir_soil(:,:)
            !alb_vis_lnd_pfts(:,:,ipft)            = alb_vis_lnd(:,:)
            !alb_nir_lnd_pfts(:,:,ipft)            = alb_nir_lnd(:,:)
            !alb_vis_can_pfts(:,:,ipft)            = alb_vis_can(:,:)
            !alb_nir_can_pfts(:,:,ipft)            = alb_nir_can(:,:)
            faPAR_pfts(:,:,ipft)                  = faPAR(:,:)
            soil_reflectivity_par_pfts(:,:,ipft)  = soil_reflectivity_par(:,:)
            ! apar_nacc_pfts(:,:,ipft)              = apar_nacc(:,:)

            faPAR_cl_pfts(:,:,:,ipft)             = faPAR_cl(:,:,:)
            apar_per_lai_cl_pfts(:,:,:,ipft)      = apar_per_lai_cl(:,:,:)
            ! apar_per_lai_cl_nacc_pfts(:,:,:,ipft) = apar_per_lai_cl_nacc(:,:,:)
            ! apar_per_lai_cl_tavg_pfts(:,:,:,ipft) = apar_per_lai_cl_tavg(:,:,:)

          END DO ! ipft

          ! Calculate the average over all lost areas
          !dummy(:,:) = SUM(sw_srf_net_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          !sw_srf_net_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(swvis_srf_net_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          !swvis_srf_net_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(swnir_srf_net_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)  / fract_lost_sum(:,:)
          !swnir_srf_net_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(lw_srf_net_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)     / fract_lost_sum(:,:)
          !lw_srf_net_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(rad_srf_net_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)    / fract_lost_sum(:,:)
          !rad_srf_net_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(alb_vis_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
          !alb_vis_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(alb_nir_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)        / fract_lost_sum(:,:)
          !alb_nir_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(alb_vis_soil_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)   / fract_lost_sum(:,:)
          !alb_vis_soil_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(alb_nir_soil_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)   / fract_lost_sum(:,:)
          !alb_nir_soil_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(alb_vis_lnd_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)    / fract_lost_sum(:,:)
          !alb_vis_lnd_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(alb_nir_lnd_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)    / fract_lost_sum(:,:)
          !alb_nir_lnd_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(alb_vis_can_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)    / fract_lost_sum(:,:)
          !alb_vis_can_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          !dummy(:,:) = SUM(alb_nir_can_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)    / fract_lost_sum(:,:)
          !alb_nir_can_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(faPAR_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)          / fract_lost_sum(:,:)
          faPAR_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(soil_reflectivity_par_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          soil_reflectivity_par_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          ! dummy(:,:) = SUM(apar_nacc_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)      / fract_lost_sum(:,:)
          ! apar_nacc_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )

          DO icanopy=1,ncanopy
            dummy(:,:) = SUM(faPAR_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)             / fract_lost_sum(:,:)
            faPAR_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            dummy(:,:) = SUM(apar_per_lai_cl_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)      / fract_lost_sum(:,:)
            apar_per_lai_cl_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            ! dummy(:,:) = SUM(apar_per_lai_cl_nacc_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
            ! apar_per_lai_cl_nacc_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
            ! dummy(:,:) = SUM(apar_per_lai_cl_tavg_pfts(:,icanopy,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
            ! apar_per_lai_cl_tavg_lost_avg(:,icanopy,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          END DO

          ! Now we can change the variables on all pfts
          CALL relocate_vars_to_new_pft_fractions( &
            & d_fract_grow_pfts(:,:,:),            &
            & f_old_to_new_grow_pfts(:,:,:),       &
            & f_d_to_new_grow_pfts(:,:,:),         &
            ! avgs
            !& sw_srf_net_lost_avg       = sw_srf_net_lost_avg(:,:,:),    &
            !& swvis_srf_net_lost_avg    = swvis_srf_net_lost_avg(:,:,:), &
            !& swnir_srf_net_lost_avg    = swnir_srf_net_lost_avg(:,:,:), &
            !& lw_srf_net_lost_avg       = lw_srf_net_lost_avg(:,:,:),    &
            !& rad_srf_net_lost_avg      = rad_srf_net_lost_avg(:,:,:),   &
            !& alb_vis_lost_avg          = alb_vis_lost_avg(:,:,:),       &
            !& alb_nir_lost_avg          = alb_nir_lost_avg(:,:,:),       &
            !& alb_vis_soil_lost_avg     = alb_vis_soil_lost_avg(:,:,:),  &
            !& alb_nir_soil_lost_avg     = alb_nir_soil_lost_avg(:,:,:),  &
            !& alb_vis_lnd_lost_avg      = alb_vis_lnd_lost_avg(:,:,:),   &
            !& alb_nir_lnd_lost_avg      = alb_nir_lnd_lost_avg(:,:,:),   &
            !& alb_vis_can_lost_avg      = alb_vis_can_lost_avg(:,:,:),   &
            !& alb_nir_can_lost_avg      = alb_nir_can_lost_avg(:,:,:),   &
            & faPAR_lost_avg            = faPAR_lost_avg(:,:,:),         &
            & soil_reflectivity_par_lost_avg     = soil_reflectivity_par_lost_avg(:,:,:), &
            ! & apar_nacc_lost_avg        = apar_nacc_lost_avg(:,:,:),     &
            ! pfts
            !& sw_srf_net    = sw_srf_net_pfts(:,:,:),    &
            !& swvis_srf_net = swvis_srf_net_pfts(:,:,:), &
            !& swnir_srf_net = swnir_srf_net_pfts(:,:,:), &
            !& lw_srf_net    = lw_srf_net_pfts(:,:,:),    &
            !& rad_srf_net   = rad_srf_net_pfts(:,:,:),   &
            !& alb_vis       = alb_vis_pfts(:,:,:),       &
            !& alb_nir       = alb_nir_pfts(:,:,:),       &
            !& alb_vis_soil  = alb_vis_soil_pfts(:,:,:),  &
            !& alb_nir_soil  = alb_nir_soil_pfts(:,:,:),  &
            !& alb_vis_lnd   = alb_vis_lnd_pfts(:,:,:),   &
            !& alb_nir_lnd   = alb_nir_lnd_pfts(:,:,:),   &
            !& alb_vis_can   = alb_vis_can_pfts(:,:,:),   &
            !& alb_nir_can   = alb_nir_can_pfts(:,:,:),   &
            & faPAR         = faPAR_pfts(:,:,:),         &
            & soil_reflectivity_par = soil_reflectivity_par_pfts(:,:,:) &
            ! & apar_nacc     = apar_nacc_pfts(:,:,:)      &
            & )

          DO icanopy=1,ncanopy
            CALL relocate_vars_to_new_pft_fractions( &
              & d_fract_grow_pfts(:,:,:),            &
              & f_old_to_new_grow_pfts(:,:,:),       &
              & f_d_to_new_grow_pfts(:,:,:),         &
              ! avgs
              & faPAR_cl_lost_avg = faPAR_cl_lost_avg(:,icanopy,:,:),                         &
              & apar_per_lai_cl_lost_avg = apar_per_lai_cl_lost_avg(:,icanopy,:,:),           &
              ! & apar_per_lai_cl_nacc_lost_avg = apar_per_lai_cl_nacc_lost_avg(:,icanopy,:,:), &
              ! & apar_per_lai_cl_tavg_lost_avg = apar_per_lai_cl_tavg_lost_avg(:,icanopy,:,:), &
              ! pfts
              & faPAR_cl = faPAR_cl_pfts(:,icanopy,:,:),                         &
              & apar_per_lai_cl = apar_per_lai_cl_pfts(:,icanopy,:,:)            &
              ! & apar_per_lai_cl_nacc = apar_per_lai_cl_nacc_pfts(:,icanopy,:,:), &
              ! & apar_per_lai_cl_tavg = apar_per_lai_cl_tavg_pfts(:,icanopy,:,:)  &
              & )
          END DO

          ! Hand down the results to child tiles
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(RAD_)%p)
                  CLASS IS (t_rad_memory)
                    rad__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_rad_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            !rad__mem%sw_srf_net%ptr(:,:)    = sw_srf_net_pfts(:,:,ipft)
            !rad__mem%swvis_srf_net%ptr(:,:) = swvis_srf_net_pfts(:,:,ipft)
            !rad__mem%swnir_srf_net%ptr(:,:) = swnir_srf_net_pfts(:,:,ipft)
            !rad__mem%lw_srf_net%ptr(:,:)    = lw_srf_net_pfts(:,:,ipft)
            !rad__mem%rad_srf_net%ptr(:,:)   = rad_srf_net_pfts(:,:,ipft)
            !rad__mem%alb_vis%ptr(:,:)       = alb_vis_pfts(:,:,ipft)
            !rad__mem%alb_nir%ptr(:,:)       = alb_nir_pfts(:,:,ipft)
            !rad__mem%alb_vis_soil%ptr(:,:)  = alb_vis_soil_pfts(:,:,ipft)
            !rad__mem%alb_nir_soil%ptr(:,:)  = alb_nir_soil_pfts(:,:,ipft)
            !rad__mem%alb_vis_lnd%ptr(:,:)   = alb_vis_lnd_pfts(:,:,ipft)
            !rad__mem%alb_nir_lnd%ptr(:,:)   = alb_nir_lnd_pfts(:,:,ipft)
            !rad__mem%alb_vis_can%ptr(:,:)   = alb_vis_can_pfts(:,:,ipft)
            !rad__mem%alb_nir_can%ptr(:,:)   = alb_nir_can_pfts(:,:,ipft)
            rad__mem%faPAR%ptr(:,:)                  = faPAR_pfts(:,:,ipft)
            rad__mem%soil_reflectivity_par%ptr(:,:)  = soil_reflectivity_par_pfts(:,:,ipft)
            ! rad__mem%apar_nacc%ptr(:,:)              = apar_nacc_pfts(:,:,ipft)
            !
            rad__mem%faPAR_cl%ptr(:,:,:)             = faPAR_cl_pfts(:,:,:,ipft)
            rad__mem%apar_per_lai_cl%ptr(:,:,:)      = apar_per_lai_cl_pfts(:,:,:,ipft)
            ! rad__mem%apar_per_lai_cl_nacc%ptr(:,:,:) = apar_per_lai_cl_nacc_pfts(:,:,:,ipft)
            ! rad__mem%apar_per_lai_cl_tavg%ptr(:,:,:) = apar_per_lai_cl_tavg_pfts(:,:,:,ipft)

          END DO ! ipft
      END IF ! RAD active


      IF (dsl4jsb_Config(TURB_)%active == .TRUE.) THEN
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(TURB_)%p)
                  CLASS IS (t_turb_memory)
                    turb__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_turb_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            fact_q_air          => turb__mem%fact_q_air%ptr(:,:)            ! inout
            fact_qsat_srf       => turb__mem%fact_qsat_srf%ptr(:,:)         ! inout
            fact_qsat_trans_srf => turb__mem%fact_qsat_trans_srf%ptr(:,:)   ! inout

            ! Collect
            fact_q_air_pfts(:,:,ipft)           = fact_q_air(:,:)
            fact_qsat_srf_pfts(:,:,ipft)        = fact_qsat_srf(:,:)
            fact_qsat_trans_srf_pfts(:,:,ipft)  = fact_qsat_trans_srf(:,:)

          END DO ! ipft

          ! Calculate the average over all lost areas
          dummy(:,:) = SUM(fact_q_air_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)          / fract_lost_sum(:,:)
          fact_q_air_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(fact_qsat_srf_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3)       / fract_lost_sum(:,:)
          fact_qsat_srf_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )
          dummy(:,:) = SUM(fact_qsat_trans_srf_pfts(:,:,:) * d_fract_shrink_pfts(:,:,:), DIM=3) / fract_lost_sum(:,:)
          fact_qsat_trans_srf_lost_avg(:,:,:)= SPREAD( dummy(:,:), DIM=2, NCOPIES=npft )

          ! Now we can change the variables on all pfts
          CALL relocate_vars_to_new_pft_fractions(   &
            & d_fract_grow_pfts(:,:,:),              &
            & f_old_to_new_grow_pfts(:,:,:),         &
            & f_d_to_new_grow_pfts(:,:,:),           &
            ! avgs
            & fact_q_air_lost_avg          = fact_q_air_lost_avg(:,:,:),          &
            & fact_qsat_srf_lost_avg       = fact_qsat_srf_lost_avg(:,:,:),       &
            & fact_qsat_trans_srf_lost_avg = fact_qsat_trans_srf_lost_avg(:,:,:), &
            ! pfts
            & fact_q_air = fact_q_air_pfts(:,:,:),                  &
            & fact_qsat_srf = fact_qsat_srf_pfts(:,:,:),            &
            & fact_qsat_trans_srf = fact_qsat_trans_srf_pfts(:,:,:) &
            & )

          ! Hand down the results to child tiles
          DO ipft=1,npft
            ! Point to the correct pft
            IF (ipft == 1) THEN
              ptr_tile => tile%Get_first_child_tile()
            ELSE
              ptr_tile => ptr_tile%Get_next_sibling_tile()
            END IF

            ! Point to the corresponding variables of this pft
            SELECT TYPE (ptr_tile)
              CLASS IS (t_jsb_tile_abstract)
              SELECT TYPE (memory => ptr_tile%mem(TURB_)%p)
                  CLASS IS (t_turb_memory)
                    turb__mem => memory
                  CLASS DEFAULT
                    CALL finish('Invalid type of ptr_tile for CLASS IS (t_turb_memory).')
                  END SELECT
              CLASS DEFAULT
                CALL finish('Invalid type of ptr_tile for CLASS IS (t_jsb_tile_abstract).')
            END SELECT

            turb__mem%fact_q_air%ptr(:,:)           = fact_q_air_pfts(:,:,ipft)
            turb__mem%fact_qsat_srf%ptr(:,:)        = fact_qsat_srf_pfts(:,:,ipft)
            turb__mem%fact_qsat_trans_srf%ptr(:,:)  = fact_qsat_trans_srf_pfts(:,:,ipft)

          END DO ! ipft
      END IF ! TURB active


      ! Remember the current pft-fractions for the next year
      fract_pfts_before_dynvegstep(:,:,:) = fract_pfts_currentyear_maps(:,:,:)


      ! Deallocate variables
      DEALLOCATE( d_fract_pfts, d_fract_grow_pfts, d_fract_shrink_pfts, fract_new_grow_pfts, f_old_to_new_grow_pfts, &
        & f_d_to_new_grow_pfts )
      DEALLOCATE( fract_lost_sum, dummy )

      IF (dsl4jsb_Config(CARBON_)%active == .TRUE.) THEN
        ! DEALLOCATE( C_pools_sum_before, C_pools_sum_after, C_pools_sum_check,   &
        DEALLOCATE( C_pools_sum_before, C_pools_sum_after,                                                                        &
          & c_green_pfts, c_reserve_pfts,  c_woods_pfts, c_crop_harvest_pfts, c_acid_ag1_pfts, c_water_ag1_pfts,                  &
          & c_ethanol_ag1_pfts, c_nonsoluble_ag1_pfts, c_acid_bg1_pfts, c_water_bg1_pfts, c_ethanol_bg1_pfts,                     &
          & c_nonsoluble_bg1_pfts, c_humus_1_pfts, c_acid_ag2_pfts, c_water_ag2_pfts, c_ethanol_ag2_pfts, c_nonsoluble_ag2_pfts,  &
          & c_acid_bg2_pfts, c_water_bg2_pfts, c_ethanol_bg2_pfts, c_nonsoluble_bg2_pfts, c_humus_2_pfts, c_bg_sum_pfts,          &
          & c_ag_sum_1_pfts, c_litter_green_ag_pfts, c_litter_green_bg_pfts, c_litter_wood_ag_pfts, c_litter_wood_bg_pfts,        &
          & NPP_pot_sum_pfts, GPP_sum_pfts, alpha_pfts, alpha_sum_pfts, alpha_yDayMean_pfts, soil_respiration_pfts,               &
          & NPP_flux_correction_pfts, NPP_pot_yDayMean_pfts, NPP_act_yDayMean_pfts, excess_NPP_pfts, GPP_yDayMean_pfts,           &
          & root_exudates_pfts, cflux_c_green_2_herb_pfts, cflux_herb_2_littergreen_pfts, cflux_herb_2_atm_pfts,                  &
          & cflux_dist_greenreserve_2_soil_pfts, cflux_dist_woods_2_soil_pfts, cflux_fire_all_2_atm_pfts, current_max_green_pfts, &
          & veg_carbon_at_max_green_pfts )
        DEALLOCATE( c_green_lost_avg, c_reserve_lost_avg, c_woods_lost_avg, c_crop_harvest_lost_avg, c_acid_ag1_lost_avg,         &
          & c_water_ag1_lost_avg, c_ethanol_ag1_lost_avg, c_nonsoluble_ag1_lost_avg, c_acid_bg1_lost_avg, c_water_bg1_lost_avg,   &
          & c_ethanol_bg1_lost_avg, c_nonsoluble_bg1_lost_avg, c_humus_1_lost_avg, c_acid_ag2_lost_avg, c_water_ag2_lost_avg,     &
          & c_ethanol_ag2_lost_avg, c_nonsoluble_ag2_lost_avg, c_acid_bg2_lost_avg, c_water_bg2_lost_avg, c_ethanol_bg2_lost_avg, &
          & c_nonsoluble_bg2_lost_avg, c_humus_2_lost_avg, c_bg_sum_lost_avg, c_ag_sum_1_lost_avg, c_litter_green_ag_lost_avg,    &
          & c_litter_green_bg_lost_avg, c_litter_wood_ag_lost_avg, c_litter_wood_bg_lost_avg,                                     &
          & NPP_pot_sum_lost_avg, GPP_sum_lost_avg, alpha_lost_avg, alpha_sum_lost_avg, alpha_yDayMean_lost_avg,                  &
          & soil_respiration_lost_avg, NPP_flux_correction_lost_avg, NPP_pot_yDayMean_lost_avg, NPP_act_yDayMean_lost_avg,        &
          & excess_NPP_lost_avg, GPP_yDayMean_lost_avg, root_exudates_lost_avg, cflux_c_green_2_herb_lost_avg,                    &
          & cflux_herb_2_littergreen_lost_avg, cflux_herb_2_atm_lost_avg,  cflux_dist_greenreserve_2_soil_lost_avg,               &
          & cflux_dist_woods_2_soil_lost_avg, cflux_fire_all_2_atm_lost_avg, current_max_green_lost_avg,                          &
          & veg_carbon_at_max_green_lost_avg )
      END IF

      IF (dsl4jsb_Config(ASSIMI_)%active == .TRUE.) THEN
        DEALLOCATE( gross_assimilation_ca_pfts, gross_assimilation_pfts, dark_respiration_ca_pfts, dark_respiration_pfts,         &
          & NPP_pot_rate_ca_pfts, canopy_cond_pfts, cover_fract_pot_pfts )
        DEALLOCATE( assimi_scaling_fact_cl_pfts, canopy_cond_cl_pfts, gross_assimilation_cl_pfts, dark_respiration_cl_pfts,       &
          & carbox_rate_max_cl_pfts, e_transport_rate_max_cl_pfts, carbox_rate_cl_pfts, e_transport_rate_cl_pfts,                 &
          & CO2_conc_leaf_cl_pfts )
        DEALLOCATE( gross_assimilation_ca_lost_avg, gross_assimilation_lost_avg, dark_respiration_ca_lost_avg,                    &
          & dark_respiration_lost_avg, NPP_pot_rate_ca_lost_avg, canopy_cond_lost_avg,            &
          & cover_fract_pot_lost_avg )
        DEALLOCATE( assimi_scaling_fact_cl_lost_avg, canopy_cond_cl_lost_avg, gross_assimilation_cl_lost_avg,                     &
          & dark_respiration_cl_lost_avg, carbox_rate_max_cl_lost_avg, e_transport_rate_max_cl_lost_avg, carbox_rate_cl_lost_avg, &
          & e_transport_rate_cl_lost_avg, CO2_conc_leaf_cl_lost_avg )
      END IF

      IF (dsl4jsb_Config(DISTURB_)%active == .TRUE.) THEN
        DEALLOCATE( cconservation_fire_pfts, cconservation_wind_pfts, burned_fract_pfts, damaged_fract_pfts )
        DEALLOCATE( cconservation_fire_lost_avg, cconservation_wind_lost_avg, burned_fract_lost_avg, damaged_fract_lost_avg )
      END IF

      IF (dsl4jsb_Config(PHENO_)%active == .TRUE.) THEN
        DEALLOCATE( fract_fpc_pfts, previous_day_NPP_pot_rate_ca_pfts, day_NPP_sum_pfts )
        DEALLOCATE( fract_fpc_lost_avg, previous_day_NPP_pot_rate_ca_lost_avg, day_NPP_sum_lost_avg )
      END IF

      IF (dsl4jsb_Config(RAD_)%active == .TRUE.) THEN
        ! DEALLOCATE( faPAR_pfts, soil_reflectivity_par_pfts, apar_nacc_pfts )
        ! DEALLOCATE( faPAR_cl_pfts, apar_per_lai_cl_pfts, apar_per_lai_cl_nacc_pfts, apar_per_lai_cl_tavg_pfts )
        ! DEALLOCATE( faPAR_lost_avg, soil_reflectivity_par_lost_avg, apar_nacc_lost_avg )
        ! DEALLOCATE( faPAR_cl_lost_avg, apar_per_lai_cl_lost_avg, apar_per_lai_cl_nacc_lost_avg,apar_per_lai_cl_tavg_lost_avg )
        DEALLOCATE( faPAR_pfts, soil_reflectivity_par_pfts)
        DEALLOCATE( faPAR_cl_pfts, apar_per_lai_cl_pfts )
        DEALLOCATE( faPAR_lost_avg, soil_reflectivity_par_lost_avg )
        DEALLOCATE( faPAR_cl_lost_avg, apar_per_lai_cl_lost_avg )
      END IF

      IF (dsl4jsb_Config(TURB_)%active == .TRUE.) THEN
        DEALLOCATE( fact_q_air_pfts, fact_qsat_srf_pfts, fact_qsat_trans_srf_pfts )
        DEALLOCATE( fact_q_air_lost_avg, fact_qsat_srf_lost_avg, fact_qsat_trans_srf_lost_avg )
      END IF

    END IF ! (is_newyear() .AND. iblk == 1)

   END SUBROUTINE dynveg_based_on_maps



  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "dynveg"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
   SUBROUTINE aggregate_dynveg(tile, options)

     USE mo_jsb_time,             ONLY: is_newyear

     CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
     TYPE(t_jsb_task_options),   INTENT(in)    :: options

     TYPE(t_jsb_model), POINTER                :: model
     CLASS(t_jsb_tile_abstract),  POINTER      :: ptr_tile

     dsl4jsb_Def_config(DYNVEG_)
     dsl4jsb_Def_memory(DYNVEG_)
     dsl4jsb_Def_config(CARBON_)
     dsl4jsb_Def_memory(CARBON_)
     dsl4jsb_Def_config(PHENO_)
     dsl4jsb_Def_memory(PHENO_)
     dsl4jsb_Def_config(ASSIMI_)
     dsl4jsb_Def_memory(ASSIMI_)
     dsl4jsb_Def_config(DISTURB_)
     dsl4jsb_Def_memory(DISTURB_)
     dsl4jsb_Def_config(RAD_)
     dsl4jsb_Def_memory(RAD_)
     dsl4jsb_Def_config(TURB_)
     dsl4jsb_Def_memory(TURB_)

     CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

     CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_dynveg'

     ! Local variables
     INTEGER  :: iblk, ics, ice

     ! Get local variables from options argument
     iblk = options%iblk
     ics  = options%ics
     ice  = options%ice

     IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

     model => Get_model(tile%owner_model_id)

     ptr_tile => tile%Get_first_child_tile()

     !dsl4jsb_Get_config(DYNVEG_)
     dsl4jsb_Get_config(CARBON_)
     dsl4jsb_Get_config(PHENO_)
     dsl4jsb_Get_config(ASSIMI_)
     dsl4jsb_Get_config(DISTURB_)
     dsl4jsb_Get_config(RAD_)
     dsl4jsb_Get_config(TURB_)

     dsl4jsb_Get_memory(DYNVEG_)

     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! Only once at first day of the year and only in the time step before the radiation time step and only once for the whole grid
    IF (is_newyear(options%current_datetime, options%dtime) .AND. &
      & is_time_ltrig_rad_m1(options%current_datetime, options%dtime, model%id)) THEN

      dsl4jsb_Aggregate_onChunk(DYNVEG_, C_pools_sum_check, weighted_by_fract)

        ! Note: Aggregate only those vars that are changed here and that are
        !       also aggregated in the processes themselfes!
        IF (dsl4jsb_Config(CARBON_)%active == .TRUE.) THEN
          dsl4jsb_Get_memory(CARBON_)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_green,           weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_reserve,         weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_woods,           weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_crop_harvest,    weighted_by_fract)
          !SIZE CLASS 1
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag1,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag1,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag1,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag1,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_bg1,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_bg1,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_bg1,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_bg1,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_1,         weighted_by_fract)
          !SIZE CLASS 2
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag2,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag2,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag2,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag2,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_bg2,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_bg2,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_bg2,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_bg2,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_2,         weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_bg_sum,          weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ag_sum_1,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_litter_green_ag, weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_litter_green_bg, weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_litter_wood_ag,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_litter_wood_bg,  weighted_by_fract)
          !
          dsl4jsb_Aggregate_onChunk(CARBON_,  NPP_pot_sum,          weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  GPP_sum,              weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  alpha,                weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  alpha_sum,            weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  alpha_yDayMean,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  soil_respiration,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  NPP_flux_correction,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  NPP_pot_yDayMean,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  NPP_act_yDayMean,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  excess_NPP,           weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  GPP_yDayMean,         weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  root_exudates,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_c_green_2_herb,           weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_herb_2_littergreen,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_herb_2_atm,               weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_greenreserve_2_soil, weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_woods_2_soil,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_fire_all_2_atm,           weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  current_max_green,              weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  veg_carbon_at_max_green,        weighted_by_fract)

          ! same for _ta variables
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_green_ta,           weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_reserve_ta,         weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_woods_ta,           weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_crop_harvest_ta,    weighted_by_fract)
          !SIZE CLASS 1
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag1_ta,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag1_ta,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag1_ta,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag1_ta,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_bg1_ta,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_bg1_ta,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_bg1_ta,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_bg1_ta,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_1_ta,         weighted_by_fract)
          !SIZE CLASS 2
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_ag2_ta,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_ag2_ta,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_ag2_ta,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_ag2_ta,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_acid_bg2_ta,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_water_bg2_ta,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_ethanol_bg2_ta,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_nonsoluble_bg2_ta,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_humus_2_ta,         weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  c_bg_sum_ta,          weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  c_ag_sum_1_ta,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_litter_green_ag_ta, weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_litter_green_bg_ta, weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_litter_wood_ag_ta,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  c_litter_wood_bg_ta,  weighted_by_fract)
          !
          !dsl4jsb_Aggregate_onChunk(CARBON_,  NPP_pot_sum_ta,          weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  GPP_sum_ta,              weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  alpha_ta,                weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  alpha_sum_ta,            weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  alpha_yDayMean_ta,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  soil_respiration_ta,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  NPP_flux_correction_ta,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  NPP_pot_yDayMean_ta,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  NPP_act_yDayMean_ta,     weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  excess_NPP_ta,           weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  GPP_yDayMean_ta,         weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  root_exudates_ta,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_c_green_2_herb_ta,           weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_herb_2_littergreen_ta,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_herb_2_atm_ta,               weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_greenreserve_2_soil_ta, weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_dist_woods_2_soil_ta,        weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(CARBON_,  cflux_fire_all_2_atm_ta,           weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  current_max_green_ta,              weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(CARBON_,  veg_carbon_at_max_green_ta,        weighted_by_fract)
        END IF ! CARBON is active


        IF (dsl4jsb_Config(ASSIMI_)%active == .TRUE.) THEN
          dsl4jsb_Get_memory(ASSIMI_)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation_ca,    weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, dark_respiration_ca,      weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, dark_respiration,         weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, NPP_pot_rate_ca,          weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond,              weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(ASSIMI_, cover_fract_pot,          weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(ASSIMI_, assimi_scaling_fact_cl,   weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_cl,           weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation_cl,    weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(ASSIMI_, dark_respiration_cl,      weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, carbox_rate_max_cl,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, e_transport_rate_max_cl,  weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, carbox_rate_cl ,          weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, e_transport_rate_cl,      weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(ASSIMI_, CO2_conc_leaf_cl,         weighted_by_fract)
        END IF ! ASSIMI is active


        IF (dsl4jsb_Config(DISTURB_)%active == .TRUE.) THEN
            dsl4jsb_Get_memory(DISTURB_)
            !dsl4jsb_Aggregate_onChunk(DISTURB_, cconservation_fire,  weighted_by_fract)
            !dsl4jsb_Aggregate_onChunk(DISTURB_, cconservation_wind,  weighted_by_fract)
            dsl4jsb_Aggregate_onChunk(DISTURB_, burned_fract,        weighted_by_fract)
            dsl4jsb_Aggregate_onChunk(DISTURB_, damaged_fract,       weighted_by_fract)
        END IF ! DISTURB is active


        IF (dsl4jsb_Config(PHENO_)%active == .TRUE.) THEN
            dsl4jsb_Get_memory(PHENO_)
            !dsl4jsb_Aggregate_onChunk(PHENO_, fract_fpc,    weighted_by_fract)
            !dsl4jsb_Aggregate_onChunk(PHENO_, previous_day_NPP_pot_rate_ca,    weighted_by_fract)
            !dsl4jsb_Aggregate_onChunk(PHENO_, cconservation_allom,    weighted_by_fract)
            !dsl4jsb_Aggregate_onChunk(PHENO_, day_NPP_sum,    weighted_by_fract)
        END IF ! PHENO is active


        IF (dsl4jsb_Config(RAD_)%active == .TRUE.) THEN
          dsl4jsb_Get_memory(RAD_)
          !dsl4jsb_Aggregate_onChunk(RAD_, sw_srf_net,     weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, swvis_srf_net,  weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, swnir_srf_net,  weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, lw_srf_net,     weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, rad_srf_net,    weighted_by_fract)

          !dsl4jsb_Aggregate_onChunk(RAD_, rad_net_lwtr,   weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, rad_net_lice,   weighted_by_fract)

          !dsl4jsb_Aggregate_onChunk(RAD_, alb_vis,        weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_nir,        weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_snow,   weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_snow,   weighted_by_fract)

          !dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_soil,   weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_soil,   weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_lnd,    weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_lnd,    weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_mineralsoil, weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_mineralsoil, weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_can,    weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_can,    weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, albedo_lwtr,    weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, albedo_lice,    weighted_by_fract)

          !dsl4jsb_Aggregate_onChunk(RAD_, faPAR,                 weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, soil_reflectivity_par, weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, apar_nacc,             weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, faPAR_cl,              weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, apar_per_lai_cl,       weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, apar_per_lai_cl_nacc,  weighted_by_fract)
          !dsl4jsb_Aggregate_onChunk(RAD_, apar_per_lai_cl_tavg,  weighted_by_fract)
        END IF ! RAD is active


        IF (dsl4jsb_Config(TURB_)%active == .TRUE.) THEN
          dsl4jsb_Get_memory(TURB_)
          !dsl4jsb_Aggregate_onChunk(TURB_, rough_h,             weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(TURB_, fact_q_air,          weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(TURB_, fact_qsat_srf,       weighted_by_fract)
          dsl4jsb_Aggregate_onChunk(TURB_, fact_qsat_trans_srf, weighted_by_fract)
        END IF ! TURB is active

     END IF ! is_newyear

     IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

   END SUBROUTINE aggregate_dynveg


#endif
END MODULE mo_dynveg_interface
