!> QUINCY agriculture process interface
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of tasks for the quincy agriculture process
!>
MODULE mo_q_agr_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_math_constants,  ONLY: eps8
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,   ONLY: t_jsb_process

  dsl4jsb_Use_processes Q_AGR_, VEG_, Q_PHENO_, Q_ASSIMI_, SB_, SPQ_, A2L_
  dsl4jsb_Use_config(Q_AGR_)
  dsl4jsb_Use_config(VEG_)
  dsl4jsb_Use_config(Q_PHENO_)
  dsl4jsb_Use_config(Q_ASSIMI_)
  dsl4jsb_Use_config(SB_)
  dsl4jsb_Use_config(SPQ_)
  dsl4jsb_Use_memory(SB_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SPQ_)
  dsl4jsb_Use_memory(Q_PHENO_)
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_memory(Q_AGR_)
  dsl4jsb_Use_memory(Q_ASSIMI_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_q_agr_tasks
  PUBLIC ::  update_cropland_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_agr_interface'

  ! ======================================================================================================= !
  !> Type definition: cropland process task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_cropland_process
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_cropland_process
    PROCEDURE, NOPASS :: Aggregate => aggregate_cropland_process
  END TYPE tsk_cropland_process

  !> Constructor interface: update_cropland_process task
  !>
  INTERFACE tsk_cropland_process
    PROCEDURE Create_task_cropland_process
  END INTERFACE tsk_cropland_process

  CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: Q_AGR_
  !>
  SUBROUTINE Register_q_agr_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_cropland_process(model_id))

  END SUBROUTINE Register_q_agr_tasks

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_cropland_process task
  !>
  FUNCTION Create_task_cropland_process(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(tsk_cropland_process::return_ptr)
    CALL return_ptr%Construct(name='cropland_process', process_id=Q_AGR_, owner_model_id=model_id)
  END FUNCTION Create_task_cropland_process

  ! ======================================================================================================= !
  !> Implementation of "update": cropland_process task
  !>
  SUBROUTINE update_cropland_process(tile, options)
    USE mo_lnd_bgcm_idx
    USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_ESTABLISHMENT_ID
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_impl_constants,      ONLY: true, false, test_false_true
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_math_constants,      ONLY: eps4, one_day
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_veg_config_class,        ONLY: get_number_of_veg_compartments
    USE mo_q_agr_process,           ONLY: calc_crop_phenology, calc_planting_flux, set_crop_n_fixation_status, &
      &                                   update_crop_growth_phase, calc_crop_allocation_factors, calc_fertiliser_application
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_lnd_bgcm_store),     POINTER       :: bgcm_store                 !< the bgcm store of this tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_lctlib_element),     POINTER       :: lctlib                     !< land-cover-type library - parameter across pft's
    INTEGER  :: lctlib_growthform, lctlib_ps_pathway
    REAL(wp) :: lctlib_gdd_req_max, lctlib_beta_soil_flush, lctlib_t_air_senescence, lctlib_beta_soil_senescence, &
      &         lctlib_cn_leaf, lctlib_np_leaf, lctlib_sla
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(Q_AGR_)
    dsl4jsb_Def_config(Q_ASSIMI_)
    dsl4jsb_Def_config(Q_PHENO_)
    dsl4jsb_Def_config(SPQ_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(VEG_)

    ! A2L_ 2D
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Real2D_onChunk :: t_air

    dsl4jsb_Def_memory(Q_AGR_)
    dsl4jsb_Real2D_onChunk :: crop_type_index
    dsl4jsb_Real2D_onChunk :: n_fertiliser_c3
    dsl4jsb_Real2D_onChunk :: n_fertiliser_c4
    dsl4jsb_Real2D_onChunk :: growth_phase
    dsl4jsb_Real2D_onChunk :: gdd_mavg
    dsl4jsb_Real2D_onChunk :: nd_crop_season
    dsl4jsb_Real2D_onChunk :: nd_crop_season_mavg

    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Real2D_onChunk :: beta_soil_gs_tphen_mavg

    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Real2D_onChunk :: growing_season
    dsl4jsb_Real2D_onChunk :: gdd
    dsl4jsb_Real2D_onChunk :: nd_dormance

    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Real2D_onChunk :: fertiliser_nh4
    dsl4jsb_Real2D_onChunk :: fertiliser_nh4_n15
    dsl4jsb_Real2D_onChunk :: fertiliser_no3
    dsl4jsb_Real2D_onChunk :: fertiliser_no3_n15
    dsl4jsb_Real2D_onChunk :: fertiliser_po4

    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    dsl4jsb_Real3D_onChunk :: soil_depth_sl

    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk :: t_air_month_mavg
    dsl4jsb_Real2D_onChunk :: t_air_tphen_mavg
    dsl4jsb_Real2D_onChunk :: t_soil_srf_tphen_mavg
    dsl4jsb_Real2D_onChunk :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk :: maint_respiration_tlabile_mavg
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real2D_onChunk :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onChunk :: leaf2sapwood_mass_ratio
    dsl4jsb_Real2D_onChunk :: leaf2root_mass_ratio
    dsl4jsb_Real2D_onChunk :: falloc_fruit_crop
    dsl4jsb_Real2D_onChunk :: active_n_fixation
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt1L2D     :: veg_establishment_mt
    dsl4jsb_Def_mt2L2D     :: veg_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                    :: ic, iblk, ics, ice, nc !< loop counter / dimensions
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb          !< Vertical grid
    INTEGER                    :: nsoil_sb               !< number of soil layers as used/defined by the SB_ process
    INTEGER                    :: nr_of_veg_bgcm_comp    !< dim for veg compartments
    REAL(wp)                   :: dtime                  !< timestep length
    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_cropland_process'

    ! ----------------------------------------------------------------------------------------------------- !
    ! Get pointers to process configs and memory
    model => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT.lctlib%CropFlag) RETURN !< do not run this routine at tiles other than croplands
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_AGR_)
    dsl4jsb_Get_config(Q_ASSIMI_)
    dsl4jsb_Get_config(Q_PHENO_)
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_config(SPQ_)
    dsl4jsb_Get_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_AGR_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    nc    = options%nc
    ics   = options%ics
    ice   = options%ice
    iblk  = options%iblk
    dtime = options%dtime
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    nr_of_veg_bgcm_comp = get_number_of_veg_compartments()
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !

    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                        ! in

    dsl4jsb_Get_memory(Q_AGR_)
    dsl4jsb_Get_var2D_onChunk(Q_AGR_, crop_type_index)            ! in
    dsl4jsb_Get_var2D_onChunk(Q_AGR_, n_fertiliser_c3)            ! in
    dsl4jsb_Get_var2D_onChunk(Q_AGR_, n_fertiliser_c4)            ! in
    dsl4jsb_Get_var2D_onChunk(Q_AGR_, growth_phase)               ! inout
    dsl4jsb_Get_var2D_onChunk(Q_AGR_, gdd_mavg)                   ! inout
    dsl4jsb_Get_var2D_onChunk(Q_AGR_, nd_crop_season)             ! inout
    dsl4jsb_Get_var2D_onChunk(Q_AGR_, nd_crop_season_mavg)        ! inout

    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tphen_mavg) ! in

    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)           ! inout
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, gdd)                      ! inout
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, nd_dormance)              ! inout

    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_var2D_onChunk(SPQ_, num_sl_above_bedrock)         ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_depth_sl)                ! in

    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_var2D_onChunk(SB_, fertiliser_nh4)                ! out
    dsl4jsb_Get_var2D_onChunk(SB_, fertiliser_nh4_n15)            ! out
    dsl4jsb_Get_var2D_onChunk(SB_, fertiliser_no3)                ! out
    dsl4jsb_Get_var2D_onChunk(SB_, fertiliser_no3_n15)            ! out
    dsl4jsb_Get_var2D_onChunk(SB_, fertiliser_po4)                ! out

    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_month_mavg)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tphen_mavg)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_soil_srf_tphen_mavg)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg) ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)    ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)    ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2sapwood_mass_ratio)      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2root_mass_ratio)         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, falloc_fruit_crop)            ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, active_n_fixation)            ! out
    ! ----------------------------------------------------------------------------------------------------- !

    lctlib_growthform = dsl4jsb_Lctlib_param(growthform)
    lctlib_ps_pathway = dsl4jsb_Lctlib_param(ps_pathway)
    lctlib_gdd_req_max = dsl4jsb_Lctlib_param(gdd_req_max)
    lctlib_beta_soil_flush = dsl4jsb_Lctlib_param(beta_soil_flush)
    lctlib_t_air_senescence = dsl4jsb_Lctlib_param(t_air_senescence)
    lctlib_beta_soil_senescence = dsl4jsb_Lctlib_param(beta_soil_senescence)
    lctlib_cn_leaf = dsl4jsb_Lctlib_param(cn_leaf)
    lctlib_np_leaf = dsl4jsb_Lctlib_param(np_leaf)
    lctlib_sla = dsl4jsb_Lctlib_param(sla)

    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)                           ! inout
    dsl4jsb_Get_mt1L2D(VEG_BGCM_ESTABLISHMENT_ID, veg_establishment_mt)         ! inout

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> - Update phenology (GDD, number of days in crop growing season), initialise growing season
    !>
    CALL calc_crop_phenology(nc,  dtime, &                                                   ! in
      &                       lctlib_beta_soil_flush, lctlib_gdd_req_max, &
      &                       lctlib_beta_soil_senescence, lctlib_t_air_senescence, &
      &                       lctlib_cn_leaf, lctlib_np_leaf, &
      &                       crop_type_index(:), t_air(:), t_air_week_mavg(:), &
      &                       t_air_month_mavg(:), beta_soil_gs_tphen_mavg(:), &
      &                       t_air_tphen_mavg(:), t_soil_srf_tphen_mavg(:), &
      &                       gpp_tlabile_mavg(:), maint_respiration_tlabile_mavg(:), &
      &                       growing_season(:), &
      &                       growth_req_n_tlabile_mavg(:), growth_req_p_tlabile_mavg(:), & ! inout
      &                       gdd(:), gdd_mavg(:), growth_phase(:), nd_dormance(:))
    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> - Determine establishment flux if this is a planting time step, and set whether the planted crop
    !>   is a N fixing crop (1) or not (0)
    !>
    CALL calc_planting_flux(nc, &                            ! in
      &                      growth_phase(:), &
      &                      veg_pool_mt(ix_seed_bed,:,:), &
      &                      veg_establishment_mt(:,:))      ! inout
    active_n_fixation(:) = set_crop_n_fixation_status(crop_type_index(:))
    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> - Update crop growth phase
    !>
    CALL update_crop_growth_phase(nc, dtime, &                                 ! in
      &                            crop_type_index(:), gdd(:), &
      &                            gdd_mavg(:), lai(:), &
      &                            nd_crop_season(:), nd_crop_season_mavg(:), & ! inout
      &                            growth_phase(:))
    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> - Calculate crop allocation factors
    !>
    CALL calc_crop_allocation_factors(nc, &
      &                                crop_type_index(:), growth_phase(:), & ! in
      &                                gdd(:), gdd_mavg(:), lai(:), &
      &                                leaf2sapwood_mass_ratio(:), &          ! out
      &                                leaf2root_mass_ratio(:), &
      &                                falloc_fruit_crop(:))
    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> - Apply fertiliser if today is a fertiliser day
    !>
    CALL calc_fertiliser_application(nc, lctlib_ps_pathway, &                      ! in
      &                               nd_crop_season(:), nd_crop_season_mavg(:), &
      &                               n_fertiliser_c3(:), n_fertiliser_c4(:), &
      &                               fertiliser_nh4(:), fertiliser_nh4_n15(:), &  ! out
      &                               fertiliser_no3(:), fertiliser_no3_n15(:), &
      &                               fertiliser_po4(:))

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_cropland_process

  ! ======================================================================================================= !
  !> Implementation of "aggregate": cropland_process task
  !>
  SUBROUTINE aggregate_cropland_process(tile, options)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_cropland_process'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_AGR_)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE aggregate_cropland_process

#endif
END MODULE mo_q_agr_interface
