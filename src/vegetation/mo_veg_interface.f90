!>
!! @brief Contains the interfaces to the vegetation processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_util,               ONLY: int2string

  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes VEG, PHENO, SOIL, SEB, ATM2LAND, RAD, HYDRO !, ASSIMI
  dsl4jsb_Use_config(VEG)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)

  dsl4jsb_Use_memory(atm2land)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(PHENO)
  dsl4jsb_Use_memory(RAD)
  dsl4jsb_Use_memory(HYDRO)
  !dsl4jsb_Use_memory(ASSIMI)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_veg_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_interface'

  !> Type definition for canopy_cond_unstressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_unstressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_unstressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_unstressed
  END TYPE tsk_canopy_cond_unstressed

  !> Constructor interface for canopy_cond_unstressed task
  INTERFACE tsk_canopy_cond_unstressed
    PROCEDURE Create_task_canopy_cond_unstressed
  END INTERFACE tsk_canopy_cond_unstressed

  !> Type definition for canopy_cond_stressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_stressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_stressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_stressed
  END TYPE tsk_canopy_cond_stressed

  !> Constructor interface for canopy_cond_stressed task
  INTERFACE tsk_canopy_cond_stressed
    PROCEDURE Create_task_canopy_cond_stressed
  END INTERFACE tsk_canopy_cond_stressed

  !> Type definition for foliage_projected_cover task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_foliage_projected_cover
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_fpc
    PROCEDURE, NOPASS :: Aggregate => aggregate_fpc
  END TYPE tsk_foliage_projected_cover

  !> Constructor interface for foliage_projected_cover task
  INTERFACE tsk_foliage_projected_cover
    PROCEDURE Create_task_foliage_projected_cover
  END INTERFACE tsk_foliage_projected_cover

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_unstressed task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "canopy_cond_unstressed"
  !!
  FUNCTION Create_task_canopy_cond_unstressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_cond_unstressed::return_ptr)

    return_ptr%name = 'canopy_cond_unstressed'
    return_ptr%process_id = VEG
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_canopy_cond_unstressed
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_stressed task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "canopy_cond_stressed"
  !!
  FUNCTION Create_task_canopy_cond_stressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_cond_stressed::return_ptr)

    return_ptr%name = 'canopy_cond_stressed'
    return_ptr%process_id = VEG
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_canopy_cond_stressed
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for foliage_projected_cover task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "foliage_projected_cover"
  !!
  FUNCTION Create_task_foliage_projected_cover(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_foliage_projected_cover::return_ptr)

    return_ptr%name = 'foliage_projected_cover'
    return_ptr%process_id = VEG
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_foliage_projected_cover
  
  SUBROUTINE Register_veg_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task(tsk_canopy_cond_unstressed(model_id))
    CALL process%Register_task(tsk_canopy_cond_stressed(model_id))
    CALL process%Register_task(tsk_foliage_projected_cover(model_id))

  END SUBROUTINE Register_veg_tasks

  SUBROUTINE update_canopy_cond_unstressed(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_veg_process, ONLY: get_canopy_conductance

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(PHENO)
    dsl4jsb_Def_memory(ATM2LAND)

    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: canopy_cond

    INTEGER :: nsoil

    INTEGER :: p_pe
    INTEGER :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_canopy_cond_unstressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(PHENO)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swpar_srf_down) ! IN
    dsl4jsb_Get_var2D_onChunk(PHENO,    lai)            ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond)    ! OUT

    ! Compute (max) canopy (stomatal) conductance with no water stress
    canopy_cond(:) = get_canopy_conductance( lai(:), swpar_srf_down(:) )

  END SUBROUTINE update_canopy_cond_unstressed

  SUBROUTINE aggregate_canopy_cond_unstressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(VEG)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_canopy_cond_unstressed'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(VEG)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(VEG,    canopy_cond,    weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_unstressed

  SUBROUTINE update_canopy_cond_stressed(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_veg_process, ONLY: get_canopy_conductance, get_water_stress_factor
    USE mo_veg_util,    ONLY: get_water_in_root_zone
    USE mo_phy_schemes, ONLY: qsat

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(VEG)
    dsl4jsb_Def_config(SOIL)

    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(PHENO)
    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_memory(HYDRO)
    !dsl4jsb_Def_memory(ASSIMI)

    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: canopy_cond
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited
    dsl4jsb_Real3D_onChunk :: soil_depth_l
    dsl4jsb_Real3D_onChunk :: root_depth_l
    !dsl4jsb_Real2D_onChunk :: w_soil_column   !< Water content in soil column
    dsl4jsb_Real3D_onChunk :: w_soil_sl            !< Water content in soil layers
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl         !< Water content of soil layers at field capacity
    dsl4jsb_Real3D_onChunk :: w_soil_pwp_sl        !< Water content of soil layers at permanent wilting point
    dsl4jsb_Real2D_onChunk :: w_soil_root         !< Water content in root zone
    dsl4jsb_Real2D_onChunk :: w_soil_root_fc      !< Water content in root zone at field capacity
    dsl4jsb_Real2D_onChunk :: w_soil_root_pwp     !< Water content in root zone at permanent wilting point
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf

    INTEGER :: nsoil

    INTEGER :: p_pe
    INTEGER :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_canopy_cond_stressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(VEG)
    dsl4jsb_Get_config(SOIL)

    nsoil = dsl4jsb_Config(SOIL)%nsoil_water

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(PHENO)
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(HYDRO)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air)               ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf)           ! in
    dsl4jsb_Get_var2D_onChunk(SEB,      t)                   ! in
    dsl4jsb_Get_var3D_onChunk(SOIL,     soil_depth_l)        ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_sl)           ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_fc_sl)        ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_pwp_sl)       ! in
    dsl4jsb_Get_var2D_onChunk(PHENO,    lai)                 ! in
    dsl4jsb_Get_var3D_onChunk(VEG,      root_depth_l)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond_limited) ! out
    dsl4jsb_Get_var3D_onChunk(VEG,      root_depth_l)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_root)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_root_fc)      ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_root_pwp)     ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO,    water_stress)        ! out

    ! Get water stress

    ! Calculate water content and field capacity in the root zone
    w_soil_root   (:)  = get_water_in_root_zone(w_soil_sl     (:,:), soil_depth_l(:,:), root_depth_l(:,:))
    w_soil_root_fc(:)  = get_water_in_root_zone(w_soil_fc_sl  (:,:), soil_depth_l(:,:), root_depth_l(:,:))
    w_soil_root_pwp(:) = get_water_in_root_zone(w_soil_pwp_sl (:,:), soil_depth_l(:,:), root_depth_l(:,:))

    ! Calculate water stress factor
    ! Temporary: later use the distributed permanent wilting point from above and change get_water_stress_factor function
    water_stress(:) = get_water_stress_factor(w_soil_root(:), w_soil_root_fc(:), w_soil_root_pwp(:))

    ! Compute (actual) canopy (stomatal) conductance under water stress.
    ! R: diese canopy_cond (stressed) entspricht der canopy_conductance_limited die in JSBACH3 in update_soil gerechnet wurde.
    !    Sie wird auf jeden Fall benötigt unabhängig davon ob der Assimilations-Prozess an ist (wird dann als Input benötigt) oder nicht.
    canopy_cond_limited(:) = get_canopy_conductance(canopy_cond(:),            & ! in, unstressed canopy conductance
                                                    water_stress(:),           & ! in, water stress factor
                                                    q_air(:) >  qsat(t(:),press_srf(:)) & ! in, atmosphere saturated?
                                                    )

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_canopy_cond_stressed

  SUBROUTINE aggregate_canopy_cond_stressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(HYDRO)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_canopy_cond_stressed'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(HYDRO)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(VEG,   canopy_cond_limited, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, w_soil_root,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, w_soil_root_fc,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, w_soil_root_pwp,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO, water_stress,        weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_stressed

  SUBROUTINE update_fpc(tile, options)

    ! ---------------------------
    ! Variables

    ! Used variables
    ! sollen jetzt in dieses File hier kommen:
    !USE mo_veg_process, ONLY: canopy_albedo, get_foliage_projected_cover
    USE mo_veg_process, ONLY: get_foliage_projected_cover

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables
    REAL(wp), PARAMETER :: &
      albedo_snow_min = 0.4_wp,  & !< Minimum possible snow albedo   ! Temporary, will depend on vegetation type
      albedo_snow_max = 0.8_wp     !< Maximum possible snow albedo   ! Temporary


    INTEGER :: p_pe

    ! Declare pointers of type process __conf and __mem
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(PHENO)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc
    dsl4jsb_Real2D_onChunk :: fract_fpc_max


    INTEGER  :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_fpc'

    ! ---------------------------
    ! Go

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

    p_pe = get_my_global_mpi_id()
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(PHENO)

    dsl4jsb_Get_var2D_onChunk(PHENO, lai)            ! in
    dsl4jsb_Get_var2D_onChunk(VEG,   fract_fpc_max)  ! in
    dsl4jsb_Get_var2D_onChunk(VEG,   fract_fpc)      ! out

    ! Update foliage projected cover
    !
    ! TBD: Use clumpiness factor instead of 1._wp
    fract_fpc(:) = get_foliage_projected_cover(fract_fpc_max(:), lai(:), SPREAD(1._wp, DIM=1, NCOPIES=SIZE(lai,1)))

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_fpc

  SUBROUTINE aggregate_fpc(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(VEG)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_fpc'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(VEG)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(VEG, fract_fpc, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_fpc

#endif
END MODULE mo_veg_interface
