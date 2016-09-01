!>
!! @brief Contains the interfaces to the soil processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish
  USE mo_util,            ONLY: One_of

  USE mo_jsb_control,     ONLY: debug_on

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: GLACIER_TYPE, VEG_TYPE, BARE_TYPE, LAND_TYPE, Contains_lct

  dsl4jsb_Use_processes SOIL_, SEB_, VEG_, PHENO_, ATM2LAND_, HYDRO_
  dsl4jsb_Use_config(SOIL_)
  dsl4jsb_Use_memory(ATM2LAND_)
  dsl4jsb_Use_memory(SOIL_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_soil_tasks

  !> Type definition for soil_temp_and_hydro task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_temp_and_hydro
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_temp_and_hydro
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_temp_and_hydro
  END TYPE tsk_soil_temp_and_hydro

  !> Constructor interface for soil_temp_and_hydro task
  INTERFACE tsk_soil_temp_and_hydro
    PROCEDURE Create_task_soil_temp_and_hydro
  END INTERFACE tsk_soil_temp_and_hydro

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_temp_and_hydro task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_temp_and_hydro"
  !!
  FUNCTION Create_task_soil_temp_and_hydro(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_temp_and_hydro::return_ptr)

    return_ptr%name = 'soil_temp_and_hydro'
    return_ptr%process_id = SOIL_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_soil_temp_and_hydro
  
  SUBROUTINE Register_soil_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task(tsk_soil_temp_and_hydro(model_id))

  END SUBROUTINE Register_soil_tasks

  ! TBD: At the moment, the hydrology for soil and canopy is handled together in the soil process (as in echam6).
  ! This will be changed later.
  SUBROUTINE update_soil_temp_and_hydro(tile, options)

    USE mo_soil_process,           ONLY: update_surf_down_soil, update_surf_down_glacier, update_soiltemp, update_surf_up
    USE mo_jsb_physical_constants, ONLY: alf

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(SOIL_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ATM2LAND_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: tte_corr
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: wind_10m
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow
    dsl4jsb_Real2D_onChunk :: w_skin_soil
    dsl4jsb_Real2D_onChunk :: w_skin
    dsl4jsb_Real2D_onChunk :: fract_water
    dsl4jsb_Real2D_onChunk :: w_snow_soil
    dsl4jsb_Real2D_onChunk :: w_snow
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: w_snow_can
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: hcap_grnd_old
    dsl4jsb_Real2D_onChunk :: evapotranspiration
    dsl4jsb_Real2D_onChunk :: trans
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: vhcap
    dsl4jsb_Real2D_onChunk :: hdiff
    dsl4jsb_Real3D_onChunk :: t_soil_acoef
    dsl4jsb_Real3D_onChunk :: t_soil_bcoef
    dsl4jsb_Real3D_onChunk :: t_soil
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: grnd_hflx_old
    dsl4jsb_Real2D_onChunk :: oro_stddev
    dsl4jsb_Real2D_onChunk :: snowmelt
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: pme_glacier
    dsl4jsb_Real2D_onChunk :: w_soil_column
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: w_soil_pwp_sl
    dsl4jsb_Real2D_onChunk :: soil_depth
    dsl4jsb_Real2D_onChunk :: w_soil_root
    dsl4jsb_Real2D_onChunk :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk :: root_depth
    dsl4jsb_Real2D_onChunk :: shc
    dsl4jsb_Real2D_onChunk :: mpot
    dsl4jsb_Real2D_onChunk :: bclapp
    dsl4jsb_Real2D_onChunk :: vpor
    dsl4jsb_Real2D_onChunk :: psdi
    dsl4jsb_Real2D_onChunk :: w_soil_overflow
    dsl4jsb_Real2D_onChunk :: w_glac


    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & t_unfilt_corrected, &
      & wsr_max, &
      & evapotrans_no_snow_skin, &
      & snow_acc, &
      & snow_melt_acc, &
      & glacier_runoff_acc, &
      & melt_water_excess, &
!vg       runoff_hd, &
!vg       drainage_hd, &
      & skin_res_evap, &
      & reduced_evap

    LOGICAL :: &
      is_glacier(options%nc)

    INTEGER :: jllog

    LOGICAL  :: lstart
    LOGICAL  :: l_tpe_open = .FALSE.
    LOGICAL  :: l_tpe_closed = .FALSE.
    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime, steplen

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_temp_and_hydro'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(SOIL_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

!!$    soil_conf => Get_soil_config(model)
    dsl4jsb_Get_config(SOIL_)

!!$    IF (model%config%usecase == 'jsbach_lite_tpe_open')   l_tpe_open = .TRUE.
!!$    IF (model%config%usecase == 'jsbach_lite_tpe_closed') l_tpe_closed = .TRUE.

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SEB_)
    IF (One_of(tile%lcts(1)%id, (/VEG_TYPE,LAND_TYPE/)) > 0) THEN
      dsl4jsb_Get_memory(VEG_)
      dsl4jsb_Get_memory(PHENO_)
    END IF
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, t_air)         ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, q_air)         ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, wind_10m)      ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, rain)          ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, snow)          ! in

    dsl4jsb_Get_var2D_onChunk(SOIL_,     hcap_grnd)     ! out
    dsl4jsb_Get_var2D_onChunk(SOIL_,     hcap_grnd_old) ! out
    dsl4jsb_Get_var2D_onChunk(SOIL_,     vhcap)         ! in
    dsl4jsb_Get_var2D_onChunk(SOIL_,     hdiff)         ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     t_soil_acoef)  ! out
    dsl4jsb_Get_var3D_onChunk(SOIL_,     t_soil_bcoef)  ! out
    dsl4jsb_Get_var3D_onChunk(SOIL_,     t_soil)        ! out
    dsl4jsb_Get_var2D_onChunk(SOIL_,     grnd_hflx)     ! out
    dsl4jsb_Get_var2D_onChunk(SOIL_,     grnd_hflx_old) ! out
    dsl4jsb_Get_var2D_onChunk(SOIL_,     soil_depth)    ! in
    dsl4jsb_Get_var2D_onChunk(SOIL_,     shc)           ! in
    dsl4jsb_Get_var2D_onChunk(SOIL_,     mpot)          ! in
    dsl4jsb_Get_var2D_onChunk(SOIL_,     bclapp)        ! in
    dsl4jsb_Get_var2D_onChunk(SOIL_,     vpor)          ! in
    dsl4jsb_Get_var2D_onChunk(SOIL_,     psdi)          ! in

    dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_skin)        ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_water)   ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_snow)        ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_snow_soil)   ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotranspiration) ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    snowmelt)      ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    runoff)        ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    drainage)      ! out

    dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt)      ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    tte_corr)      ! out
    dsl4jsb_Get_var2D_onChunk(SOIL_,     oro_stddev)    ! in

    IF (tile%lcts(1)%id == GLACIER_TYPE) THEN
      is_glacier(:) = .TRUE.
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_glac)        ! inout
    ELSE
      is_glacier(:) = .FALSE.
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    w_soil_sl)     ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    w_soil_fc_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    w_soil_pwp_sl) ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_soil_column) ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_skin_soil)   ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_soil_overflow) ! inout
      IF (One_of(tile%lcts(1)%id, (/VEG_TYPE ,LAND_TYPE/)) > 0) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_snow_can)    ! inout
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_soil_root)   ! inout
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    w_soil_root_fc) ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    trans)         ! in
        dsl4jsb_Get_var2D_onChunk(VEG_,   fract_fpc_max)    ! in
        dsl4jsb_Get_var2D_onChunk(VEG_,   root_depth)       ! in
      END IF
      IF (One_of(tile%lcts(1)%id, (/VEG_TYPE,LAND_TYPE/)) > 0) THEN
        dsl4jsb_Get_var2D_onChunk(PHENO_, lai)              ! in
      END IF
    END IF

    IF (ALL(w_snow(:) < 0._wp)) THEN          ! Start of experiment
      lstart = .TRUE.
      w_snow(:) = 0._wp
    ELSE
      lstart = .FALSE.
    END IF

    ! Maximum capacity of skin reservoir (soil + canopy)
    IF (One_of(tile%lcts(1)%id, (/VEG_TYPE,LAND_TYPE/)) > 0) THEN
      wsr_max(:) = dsl4jsb_Config(SOIL_)%wsr_capacity * (1._wp + lai(:) * fract_fpc_max(:))
    ELSE
      wsr_max(:) = dsl4jsb_Config(SOIL_)%wsr_capacity
    END IF

    !IF (lstart) THEN
      evapotrans_no_snow_skin(:) = 0._wp
      snow_acc(:) = 0._wp
      snow_melt_acc(:) = 0._wp
      glacier_runoff_acc(:) = 0._wp
      melt_water_excess(:) = 0._wp
    !END IF

    t_unfilt_corrected(:) = t_unfilt(:)

    IF (tile%lcts(1)%id == GLACIER_TYPE) THEN
      CALL update_surf_down_glacier( &
        & dtime,                      & !< in
        & lstart,                     & !< in
        & t_unfilt(:),                & !< in
        & fract_snow(:),              & !< in
        & hcap_grnd(:),               & !< in
        & evapotranspiration(:),      & !< in
        & evapopot(:),                & !< in
        & rain(:),                    & !< in
        & snow(:),                    & !< in
        & w_glac(:),                  & !< inout
        & evapotrans_no_snow_skin(:), & !< out
        & runoff(:),                  & !< out
        !& snow_acc(:),                &
        !& snow_melt_acc(:),           &
        !& glacier_runoff_acc(:),      &
        & snowmelt(:),                & !< out
        & tte_corr(:)                 & !< out
        & )
      drainage(:) = 0._wp
      w_snow_soil(:) = 0._wp
    ELSE
      CALL update_surf_down_soil(     &
        & dtime,                      & !< in
        & lstart,                     & !< in
        & q_air(:),                   & !< in
        & t_unfilt(:),                & !< in
        & wind_10m(:),                & !< in
        & t_air(:),                   & !< in
        & w_skin_soil(:),             & !< inout
        & fract_water(:),             & !< in
        & wsr_max(:),                 & !< in
        & w_snow_soil(:),             & !< inout
        & fract_snow(:),              & !< in
        & w_snow_can(:),              & !< inout
        & hcap_grnd(:),               & !< in
        & evapotranspiration(:),      & !< in
        & evapopot(:),                & !< in
        & evapotrans_no_snow_skin(:), & !< out
        & rain(:),                    & !< in
        & snow(:),                    & !< in
        !& snow_acc(:),                &
        !& snow_melt_acc(:),           &
        !& glacier_runoff_acc(:),      &
        & melt_water_excess(:),       & !< out
        & snowmelt(:),                & !< out
        & tte_corr(:)                 & !< out
        & )
    END IF
    ! Save old ground heat flux and heat capacity before updating
    ! The old flux and capacity (which have been used in the implicit calculation of the new surface temperature) are
    ! given back to the atmosphere so that the surface energy balance can be computed there
    grnd_hflx_old(:) = grnd_hflx(:)
    hcap_grnd_old(:) = hcap_grnd(:)

    ! Correction of surface temperature for snow and ice melt
    IF (.NOT. lstart) THEN
      t_unfilt_corrected(:) = t_unfilt(:) - snowmelt(:) * dtime * alf / hcap_grnd(:) 
    ELSE
      t_unfilt_corrected(:) = t_unfilt(:)
    END IF

    CALL update_soiltemp( &
      & nc,                &  !< in
      & dsl4jsb_Config(SOIL_)%nsoil_water, & !< in
      & dtime,             & !< in
      & lstart,            & !< in
      & t_unfilt_corrected(:), & !< in
      & w_snow_soil(:),    & !< in
      & hdiff(:),          & !< in
      & vhcap(:),          & !< in
      & t_soil_acoef(:,:), & !< inout
      & t_soil_bcoef(:,:), & !< inout
      & t_soil(:,:),       & !< inout
      & hcap_grnd(:),      & !< out
      & grnd_hflx(:),      & !< out
      & is_glacier(:)      & !< in
      & )

    IF (tile%lcts(1)%id == GLACIER_TYPE) THEN
      ! Don't change evapotrans_no_snow_skin
      ! Runoff and drainage already computed above
      skin_res_evap(:) = 0._wp
    ELSE
      !vg:  runoff_hd and drainage_hd no longer needed
      !vg:  accumulation does not make sense (no reset)
      CALL update_surf_up( &
        & nc,                               & !< in
        & dtime,                            & !< in
        & dsl4jsb_Config(SOIL_)%nsoil_water, & !< in
        & w_skin_soil(:),                   & !< inout
        & fract_water(:),                   & !< in
        & wsr_max(:),                       & !< in
        & w_soil_root(:),                   & !< inout
        & w_soil_root_fc(:),                & !< in
        & t_soil(:,1),                      & !< in
        & fract_snow(:),                    & !< in
        & oro_stddev(:),                    & !< in
        & evapotranspiration(:),            & !< in
        & evapopot(:),                      & !< in
        & evapotrans_no_snow_skin(:),       & !< inout ; TBD: diagnostic, currently not used
        & rain(:),                          & !< in
        & runoff(:),                        & !< inout
        & drainage(:),                      &  !< inout
        !vg runoff_hd(:), drainage_hd(:), &
        & melt_water_excess(:),             & !< in
        & skin_res_evap(:),                 & !< out ; TBD: diagnostic, currently not used
        & dsl4jsb_Config(SOIL_)%dz_water(1:dsl4jsb_Config(SOIL_)%nsoil_water), & !< in
        & w_soil_sl(:,:),                   & !< inout
        & w_soil_fc_sl(:,:),                & !< in
        & w_soil_pwp_sl(:,:),               & !< in
        & root_depth(:),                    & !< in
        & soil_depth(:),                    & !< in
        & shc(:),                           & !< in
        & mpot(:),                          & !< in
        & bclapp(:),                        & !< in
        & vpor(:),                          & !< in
        & trans(:),                         & !< in
        & reduced_evap(:),                  & !< out
  !!$      & SPREAD(0.318517_wp, NCOPIES=nc, DIM=1), & !Temporary
  !!$      & SPREAD(0.185245_wp, NCOPIES=nc, DIM=1), & !Temporary
        & psdi(:),                          & !< in
        & w_soil_overflow(:),               & !< inout
        & jllog,                            & !< out
        & l_tpe_open,                       & !< in
        & l_tpe_closed                      & !< in
        & )
    END IF

    SELECT CASE (tile%lcts(1)%id)
    CASE (GLACIER_TYPE)
      w_snow(:) = 0._wp
    CASE (BARE_TYPE)
      w_soil_column(:) = SUM(w_soil_sl(:,:), DIM=2)
      w_snow(:) = w_snow_soil(:)
      w_skin(:) = w_skin_soil(:)
    CASE (VEG_TYPE)
      w_soil_column(:) = SUM(w_soil_sl(:,:), DIM=2)
      w_snow(:) = w_snow_soil(:) + w_snow_can(:)
      w_skin(:) = w_skin_soil(:) ! + w_skin_can(:)
    END SELECT

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_temp_and_hydro

  SUBROUTINE aggregate_soil_temp_and_hydro(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(SOIL_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_soil_temp_and_hydro'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract", copy=.TRUE.)

    dsl4jsb_Aggregate_onChunk(SOIL_,  hcap_grnd,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SOIL_,  hcap_grnd_old,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SOIL_,  t_soil_acoef,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SOIL_,  t_soil_bcoef,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SOIL_,  t_soil,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SOIL_,  grnd_hflx,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SOIL_,  grnd_hflx_old,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_skin,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_snow,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snowmelt,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, tte_corr,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_glac,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_soil_sl,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_soil_column,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_skin_soil,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_snow_soil,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_soil_overflow, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_snow_can,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, w_soil_root,     weighted_by_fract)

    CALL weighted_by_fract%Set_option('l_all_tiles', .TRUE.)

    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,        weighted_by_fract)

    DEALLOCATE(weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_temp_and_hydro

  SUBROUTINE update_soil_2(tile, options)

    ! ---------------------------
    ! Variables

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables


    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & fract          (:)

    ! Locally allocated vectors
    !
!!$    REAL(wp) :: &
    INTEGER :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_2'

    ! ---------------------------
    ! Go

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(SOIL_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    CALL tile%Get_fraction(ics, ice, iblk, fract=fract(:))
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN

  END SUBROUTINE update_soil_2

    SUBROUTINE aggregate_soil(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_soil'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    model => Get_model(tile%owner_model_id)

  END SUBROUTINE aggregate_soil

#endif
END MODULE mo_soil_interface
