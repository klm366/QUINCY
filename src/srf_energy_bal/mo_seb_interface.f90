  !>
  !! @brief Contains the interfaces for the surface energy balance process.
  !!
  !! @author
  !!  Reiner Schnur, MPI-M Hamburg
  !!
  !! @par Revision History
  !! First version in mo_srf_interface.f90 by Reiner Schnur (2013-04-04)
  !! Separated out the surface energy balance part by Reiner Schnur (2016-03-22)
  !!
MODULE mo_seb_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_lct_class,      ONLY: LAND_TYPE, LAKE_TYPE, Contains_lct
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes SEB, HYDRO, SOIL
  dsl4jsb_Use_config(SEB)

  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(HYDRO)
  dsl4jsb_Use_memory(SOIL)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_seb_tasks !,t_seb_process

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_energy
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_energy
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_energy
  END TYPE tsk_surface_energy

  INTERFACE tsk_surface_energy
    PROCEDURE Create_task_surface_energy
  END INTERFACE tsk_surface_energy

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snowmelt_correction
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snowmelt_correction
    PROCEDURE, NOPASS :: Aggregate => aggregate_snowmelt_correction
  END TYPE tsk_snowmelt_correction

  INTERFACE tsk_snowmelt_correction
    PROCEDURE Create_task_snowmelt_correction
  END INTERFACE tsk_snowmelt_correction

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_asselin_filter
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_asselin
    PROCEDURE, NOPASS :: Aggregate => aggregate_asselin
  END TYPE tsk_asselin_filter

  INTERFACE tsk_asselin_filter
    PROCEDURE Create_task_asselin_filter
  END INTERFACE tsk_asselin_filter

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_fluxes
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_fluxes
  END TYPE tsk_surface_fluxes

  INTERFACE tsk_surface_fluxes
    PROCEDURE Create_task_surface_fluxes
  END INTERFACE tsk_surface_fluxes

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_energy task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_energy"
  !!
  FUNCTION Create_task_surface_energy(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_energy::return_ptr)

    return_ptr%name = 'surface_energy'
    return_ptr%process_id = SEB
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_surface_energy

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snowmelt_correction task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "snowmelt_correction"
  !!
  FUNCTION Create_task_snowmelt_correction(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snowmelt_correction::return_ptr)

    return_ptr%name = 'snowmelt_correction'
    return_ptr%process_id = SEB
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_snowmelt_correction

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for asselin_filter task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "asselin_filter"
  !!
  FUNCTION Create_task_asselin_filter(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_asselin_filter::return_ptr)

    return_ptr%name = 'asselin_filter'
    return_ptr%process_id = SEB
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_asselin_filter

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_energy fluxes
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_fluxes"
  !!
  FUNCTION Create_task_surface_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_fluxes::return_ptr)

    return_ptr%name = 'surface_fluxes'
    return_ptr%process_id = SEB
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_surface_fluxes
  
  ! =======================================================================================================
  !> Register tasks for surface energy process
  !!
  !! @param[in,out] this      Instance of surface energy process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_seb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_surface_energy(model_id))
    CALL this%Register_task(tsk_snowmelt_correction(model_id))
    CALL this%Register_task(tsk_asselin_filter(model_id))
    CALL this%Register_task(tsk_surface_fluxes(model_id))

  END SUBROUTINE Register_seb_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "surface energy"
  !! Task "update_surface_energy" calculates the new surface temperature from the surface energy balance (sensible and latent heat, 
  !! net radiation, ground heat).
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_energy(tile, options)

    USE mo_seb_land, ONLY: update_surface_energy_land
    USE mo_seb_lake, ONLY: update_surface_energy_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk, ics, ice

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy'

    IF (.NOT. tile%Is_process_active(SEB)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      CALL update_surface_energy_lake(tile, options)
    ELSE IF (Contains_lct(tile%lcts, LAND_TYPE)) THEN
      CALL update_surface_energy_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_surface_energy

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_energy"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_surface_energy(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_surface_energy'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SEB)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SEB, t_unfilt,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, t_eff,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, qsat_star,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, t_unfilt,      weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SEB, lw_srf_net,     weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SEB, sw_srf_net,     weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SEB, swvis_srf_net,  weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SEB, swnir_srf_net,  weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(SEB, t_lwtr,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, t_lice,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, s_wtr,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, s_ice,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, qsat_lwtr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, qsat_lice, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_energy

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "snowmelt_correction"
  !! Task "snowmelt_correction" applies a correction for t_unfilt due to snow and ice melting
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snowmelt_correction(tile, options)

    USE mo_jsb_physical_constants, ONLY: alf

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(HYDRO)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: snowmelt
    dsl4jsb_Real2D_onChunk :: hcap_grnd_old

    INTEGER  :: iblk, ics, ice
    REAL(wp) :: dtime

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_snowmelt_correction'

    IF (.NOT. tile%Is_process_active(SEB)) RETURN

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    dtime = options%dtime

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(HYDRO)
    IF (Contains_lct(tile%lcts, LAND_TYPE)) THEN
      dsl4jsb_Get_memory(SOIL)
    END IF

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      ! No correction for lake tile
    ELSE IF (Contains_lct(tile%lcts, LAND_TYPE)) THEN
      dsl4jsb_Get_var2D_onChunk(SEB,      t_unfilt)      ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO,    snowmelt)      ! in
      dsl4jsb_Get_var2D_onChunk(SOIL,     hcap_grnd_old) ! in

      IF (ALL(hcap_grnd_old(:) > EPSILON(1._wp))) &  ! .NOT. lstart
        & t_unfilt(:) = t_unfilt(:) - snowmelt(:) * dtime * alf / hcap_grnd_old(:) 
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_snowmelt_correction

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snowmelt_correction"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_snowmelt_correction(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_snowmelt_correction'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SEB)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SEB, t_unfilt, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snowmelt_correction
  
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "asselin_filter"
  !! Task "asselin_filter" calculates applies the Asselin time filter to the new surface temperature, if applicable
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_asselin(tile, options)

    USE mo_seb_land, ONLY: update_asselin_land

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk, ics, ice

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_asselin'

    IF (.NOT. tile%Is_process_active(SEB)) RETURN

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    model => Get_model(tile%owner_model_id)

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      ! No Asselin filter for lake tile
    ELSE IF (Contains_lct(tile%lcts, LAND_TYPE)) THEN
      CALL update_asselin_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_asselin

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "asselin_filter"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_asselin(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_asselin'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SEB)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SEB, t,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB, t_filt, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_asselin
  
  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "surface fluxes"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_fluxes(tile, options)

    USE mo_seb_land, ONLY: update_surface_fluxes_land
    USE mo_seb_lake, ONLY: update_surface_fluxes_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_fluxes'

    IF (.NOT. tile%Is_process_active(SEB)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      CALL update_surface_fluxes_lake(tile, options)
    ELSE IF (Contains_lct(tile%lcts, LAND_TYPE)) THEN
      CALL update_surface_fluxes_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_surface_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_fluxes"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.  
  SUBROUTINE aggregate_surface_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(HYDRO)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_surface_fluxes'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(HYDRO)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO, evapotranspiration, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, evapopot,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, trans,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB,   latent_hflx,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB,   sensible_hflx,      weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(HYDRO, evapo_lnd,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, evapo_wtr,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB,   latent_hflx_wtr,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB,   sensible_hflx_wtr,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, evapo_ice,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB,   latent_hflx_ice,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB,   sensible_hflx_ice,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_fluxes
  !
  ! ================================================================================================================================
  !
#endif
END MODULE mo_seb_interface
