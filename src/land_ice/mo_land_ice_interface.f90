!> Contains the interfaces to the land_ice process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Computes surface mass balance and collects fields for coupling to an external ice-sheet model (e.g. Elmer/Ice)
!>
MODULE mo_land_ice_interface
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_LAND_ICE__)

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes LAND_ICE_

  ! Use of process configurations
  dsl4jsb_Use_config(LAND_ICE_)

  ! Use of process memories
  dsl4jsb_Use_memory(LAND_ICE_)

  ! Use of other process memories for input
  dsl4jsb_Use_processes SEB_, HYDRO_
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_land_ice_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_land_ice_interface'

  !> Type definition for land_ice task: update surface mass balance
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_land_ice_update_smb
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_land_ice_smb
    PROCEDURE, NOPASS :: Aggregate => aggregate_land_ice_smb
  END TYPE tsk_land_ice_update_smb

  !> Constructor interface for land_ice smb task
  INTERFACE tsk_land_ice_update_smb
    PROCEDURE Create_task_land_ice_update_smb
  END INTERFACE tsk_land_ice_update_smb

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for land_ice_update_smb task
  !!
  FUNCTION Create_task_land_ice_update_smb(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_land_ice_update_smb::return_ptr)
    CALL return_ptr%Construct(name='land_ice_update_smb', process_id=LAND_ICE_, owner_model_id=model_id)

  END FUNCTION Create_task_land_ice_update_smb

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for land_ice process
  !!
  SUBROUTINE Register_land_ice_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task(tsk_land_ice_update_smb(model_id))

  END SUBROUTINE Register_land_ice_tasks

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "land_ice_update_smb"
  !! Collects surface mass balance components from SEB and HYDRO processes on the glacier tile
  !! and writes them into the LAND_ICE memory for later coupling to the ice-sheet model.
  !!
  SUBROUTINE update_land_ice_smb(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_land_ice_smb'

    ! Declare process memories
    dsl4jsb_Def_config(LAND_ICE_)
    dsl4jsb_Def_memory(LAND_ICE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Declare pointers to variables in memory
    ! -- LAND_ICE output
    dsl4jsb_Real2D_onChunk :: surface_temp
    dsl4jsb_Real2D_onChunk :: surface_mass_balance
    dsl4jsb_Real2D_onChunk :: ground_heat_flux
    dsl4jsb_Real2D_onChunk :: rainfall_rate
    dsl4jsb_Real2D_onChunk :: snowfall_rate

    ! -- SEB input
    dsl4jsb_Real2D_onChunk :: t_srf_unfilt
    dsl4jsb_Real2D_onChunk :: grnd_hflx

    ! -- HYDRO input
    dsl4jsb_Real2D_onChunk :: evapotrans
    dsl4jsb_Real2D_onChunk :: water_snow
    dsl4jsb_Real2D_onChunk :: precip_rain
    dsl4jsb_Real2D_onChunk :: precip_snow

    INTEGER :: ic

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(LAND_ICE_)
    IF (.NOT. dsl4jsb_Config(LAND_ICE_)%active) RETURN

    ! Get process memories
    dsl4jsb_Get_memory(LAND_ICE_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

!     ! Get input variables from SEB
!     dsl4jsb_Get_var2D_onChunk(SEB_,   t_srf_unfilt)        ! IN
!     dsl4jsb_Get_var2D_onChunk(SEB_,   grnd_hflx)           ! IN
!
!     ! Get input variables from HYDRO
!     dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)           ! IN
!     dsl4jsb_Get_var2D_onChunk(HYDRO_, precip_rain)          ! IN
!     dsl4jsb_Get_var2D_onChunk(HYDRO_, precip_snow)          ! IN
!
!     ! Get output variables in LAND_ICE
!     dsl4jsb_Get_var2D_onChunk(LAND_ICE_, surface_temp)           ! OUT
!     dsl4jsb_Get_var2D_onChunk(LAND_ICE_, surface_mass_balance)   ! OUT
!     dsl4jsb_Get_var2D_onChunk(LAND_ICE_, ground_heat_flux)       ! OUT
!     dsl4jsb_Get_var2D_onChunk(LAND_ICE_, rainfall_rate)          ! OUT
!     dsl4jsb_Get_var2D_onChunk(LAND_ICE_, snowfall_rate)          ! OUT
!
!     !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
!     DO ic = 1, nc
!       surface_temp(ic)           = t_srf_unfilt(ic)
!       ground_heat_flux(ic)       = grnd_hflx(ic)
!       rainfall_rate(ic)          = precip_rain(ic)
!       snowfall_rate(ic)          = precip_snow(ic)
!       ! SMB = snowfall - sublimation/evaporation (evapotrans is negative for sublimation)
!       ! Note: a more refined SMB calculation may be needed depending on the ice-sheet model requirements
!       surface_mass_balance(ic)   = precip_snow(ic) + evapotrans(ic)
!     END DO
!     !$ACC END PARALLEL LOOP

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_land_ice_smb

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "land_ice_update_smb"
  !!
  SUBROUTINE aggregate_land_ice_smb(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_land_ice_smb'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

!     dsl4jsb_Aggregate_onChunk(LAND_ICE_, surface_temp,         weighted_by_fract)
!     dsl4jsb_Aggregate_onChunk(LAND_ICE_, surface_mass_balance, weighted_by_fract)
!     dsl4jsb_Aggregate_onChunk(LAND_ICE_, ground_heat_flux,     weighted_by_fract)
!     dsl4jsb_Aggregate_onChunk(LAND_ICE_, rainfall_rate,        weighted_by_fract)
!     dsl4jsb_Aggregate_onChunk(LAND_ICE_, snowfall_rate,        weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_land_ice_smb

#endif
END MODULE mo_land_ice_interface
