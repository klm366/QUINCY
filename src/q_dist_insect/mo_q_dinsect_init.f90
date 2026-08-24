!> Initialize insect disturbance variables (QUINCY)
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### initialization of insect disturbance memory variables
!>
MODULE mo_q_dinsect_init
#ifndef __NO_JSBACH__

    USE mo_kind,                    ONLY: wp
    USE mo_exception,               ONLY: message, finish, message_text
    USE mo_jsb_control,             ONLY: debug_on
    USE mo_jsb_math_constants,      ONLY: eps8

    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_grid
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract

    USE mo_lnd_bgcm_idx
    USE mo_veg_config_class,        ONLY: get_number_of_veg_compartments
    USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store, get_bgcm_idx
    USE mo_lnd_bgcm_store_class,    ONLY: DINSECT_BGCM_DEADSAPWOOD_ID,   &
     &                                    DINSECT_BGCM_DEADHEARTWOOD_ID, &
     &                                    DINSECT_BGCM_DEADLEAF_ID,      &
     &                                    DINSECT_BGCM_LARVA_ID

    USE mo_jsb_process_class,     ONLY: Q_DINSECT_

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: dinsect_init

    CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_init'

  CONTAINS

    !> Initialize insect disturbance process (after memory has been set up)
    !>
    !> Only for QLAND model
    !>
  SUBROUTINE dinsect_init(tile)
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_quincy_model_config,             ONLY: QLAND
    USE mo_jsb_task_class,                  ONLY: t_jsb_task_options
    USE mo_jsb_model_class,                 ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_jsb_process_class,               ONLY: Q_DINSECT_
    USE mo_jsb_grid_class,                  ONLY: t_jsb_grid
    USE mo_jsb_grid,                        ONLY: Get_grid

    USE mo_q_dinsect_constants,               ONLY: !

    ! ----------------------------------------------------------------------------------------------------- !

    dsl4jsb_Use_config(Q_DINSECT_)
    dsl4jsb_Use_memory(Q_DINSECT_)

    IMPLICIT NONE
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile          !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),              POINTER :: model               !< the model
    TYPE(t_lnd_bgcm_store),         POINTER :: bgcm_store          !< the bgcm store of this tile
    TYPE(t_jsb_grid),               POINTER :: hgrid               !< Horizontal grid
    TYPE(t_lctlib_element),         POINTER :: lctlib              !< land-cover-type library - parameter across pft's

    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: hlp1                !< dummy variable
    INTEGER                                 :: iblk, ic            !< loop indices
    INTEGER                                 :: nproma, nblks       !< dimensions
    INTEGER                     :: dinsect_swood_idx               !< index of the dead sapwood within the bgcm store
    INTEGER                     :: dinsect_hwood_idx               !< index of the dead heartwood within the bgcm store
    INTEGER                     :: dinsect_leaf_idx                !< index of the dead leaf within the bgcm store
    INTEGER                     :: dinsect_larva_idx               !< index of the larva within the bgcm store
    REAL(wp),           POINTER :: dinsect_swood_mt_domain(:,:,:)  !< dim: elements, nc, nblks
    REAL(wp),           POINTER :: dinsect_hwood_mt_domain(:,:,:)  !< dim: elements, nc, nblks
    REAL(wp),           POINTER :: dinsect_leaf_mt_domain(:,:,:)   !< dim: elements, nc, nblks
    REAL(wp),           POINTER :: dinsect_larva_mt_domain(:,:,:)  !< dim: elements, nc, nblks
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':dinsect_init'
    ! ----------------------------------------------------------------------------------------------------- !

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DINSECT_)
    dsl4jsb_Def_memory(Q_DINSECT_)
    ! ----------------------------------------------------------------------------------------------------- !

    IF (.NOT. tile%Is_process_calculated(Q_DINSECT_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    IF (tile%lcts(1)%lib_id == 0) RETURN                           !< run this init only if the present tile is a pft
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    hgrid  => Get_grid(model%grid_id)
    nblks  =  hgrid%nblks
    nproma =  hgrid%nproma
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DINSECT_)
    dsl4jsb_Get_memory(Q_DINSECT_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(hlp1(nproma, nblks))
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    !> init '?_pool_mt' with zero below
    dinsect_swood_idx = get_bgcm_idx(bgcm_store, DINSECT_BGCM_DEADSAPWOOD_ID, tile%name, routine)
    dinsect_swood_mt_domain => bgcm_store%store_1l_2d_bgcms(bgcm_store%idx_in_store(dinsect_swood_idx))%mt_1l_2d_bgcm(:,:,:)
    dinsect_hwood_idx = get_bgcm_idx(bgcm_store, DINSECT_BGCM_DEADHEARTWOOD_ID, tile%name, routine)
    dinsect_hwood_mt_domain => bgcm_store%store_1l_2d_bgcms(bgcm_store%idx_in_store(dinsect_hwood_idx))%mt_1l_2d_bgcm(:,:,:)
    dinsect_leaf_idx = get_bgcm_idx(bgcm_store, DINSECT_BGCM_DEADLEAF_ID, tile%name, routine)
    dinsect_leaf_mt_domain => bgcm_store%store_1l_2d_bgcms(bgcm_store%idx_in_store(dinsect_leaf_idx))%mt_1l_2d_bgcm(:,:,:)
    dinsect_larva_idx = get_bgcm_idx(bgcm_store, DINSECT_BGCM_LARVA_ID, tile%name, routine)
    dinsect_larva_mt_domain => bgcm_store%store_1l_2d_bgcms(bgcm_store%idx_in_store(dinsect_larva_idx))%mt_1l_2d_bgcm(:,:,:)
    ! ----------------------------------------------------------------------------------------------------- !
    !> init dinsect_pool_mt with zero
    dinsect_swood_mt_domain(:,:,:) = 0.0_wp
    dinsect_hwood_mt_domain(:,:,:) = 0.0_wp
    dinsect_leaf_mt_domain(:,:,:) = 0.0_wp
    dinsect_larva_mt_domain(:,:,:) = 0.0_wp

    DEALLOCATE(hlp1)

  END SUBROUTINE dinsect_init

  ! ====================================================================================================== !
  !
  !> Reading of this year to be consumed by bark beetle in grid-cell area from input file
  !
  SUBROUTINE read_bb_dist_fraction(model_id, current_datetime)
    USE mo_io_units,              ONLY: filename_max
    USE mo_jsb_time_iface,        ONLY: t_datetime
    USE mo_jsb_time,              ONLY: is_time_experiment_start, get_year
    USE mo_jsb_io_netcdf,         ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_class,             ONLY: Get_model
    dsl4jsb_Use_config(Q_DINSECT_)
    dsl4jsb_Use_memory(Q_DINSECT_)
    ! USE mo_q_syl_config_class,    ONLY: CONST_HARV, REF_HARV, TRANS_HARV

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),            POINTER :: model
    TYPE(t_jsb_grid),             POINTER :: hgrid
    CLASS(t_jsb_tile_abstract),   POINTER :: tile
    CHARACTER(len=:), ALLOCATABLE         :: fullpath

    dsl4jsb_Def_config(Q_DINSECT_)
    dsl4jsb_Def_memory(Q_DINSECT_)

    dsl4jsb_Real2D_onDomain :: dinsect_bb_dist_fraction_read
    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp),           POINTER :: ptr_2D(:,:)  ! tmp pointer
    INTEGER                     :: current_year
    TYPE(t_input_file)          :: input_file
    CHARACTER(len=filename_max) :: fname
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_harvest_fraction'
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)
    CALL model%Get_top_tile(tile)

    dsl4jsb_Get_config(Q_DINSECT_)
    dsl4jsb_Get_memory(Q_DINSECT_)
    dsl4jsb_Get_var2D_onDomain(Q_DINSECT_, dinsect_bb_dist_fraction_read)

    IF (.NOT. dsl4jsb_Config(Q_DINSECT_)%dinsect_flag_bb_from_namelist) THEN
      ! No reading necessary if constant harvest
      RETURN
    END IF
    ! IF (dsl4jsb_Config(Q_DINSECT_)%harvest_scheme == CONST_HARV) THEN
    !   ! No reading necessary if constant harvest
    !   RETURN
    ! END IF

    ! IF (dsl4jsb_Config(Q_DINSECT_)%harvest_scheme == REF_HARV &
    !     & .AND. .NOT. is_time_experiment_start(current_datetime)) THEN
    !   ! In case of a reference harvest file, reading is only necessary at the beginning of the experiment
    !   RETURN
    ! END IF

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    ! IF (dsl4jsb_Config(Q_DINSECT_)%harvest_scheme == REF_HARV) THEN
    !   ! search a file not ending on a year ...
    !   WRITE (fname,'(a,a)') &
    !     & TRIM(dsl4jsb_Config(Q_DINSECT_)%harvest_filename_prefix), ".nc"
    ! ELSE
      ! search a file ending on the current year
      current_year  = get_year(current_datetime)

      !>
      !> Assertion: routine currently expects filenames with 4 digits
      !>
      IF ( current_year > 9999) THEN
        CALL finish(TRIM(routine), 'Violation of assertion: this routine currently expects filenames with 4 digits.')
      ELSE IF (current_year < 1000) THEN
        CALL message(TRIM(routine), 'Warning: this routine currently expects filenames with 4 digits, '&
                                    //' for years < 1000 leading zeros are expected.')
      END IF
      WRITE (fname,'(a,a,I4.4,a)') &
        & TRIM(dsl4jsb_Config(Q_DINSECT_)%dinsect_bb_filename_prefix), '_', current_year, ".nc"
    ! END IF
    fullpath = TRIM(dsl4jsb_Config(Q_DINSECT_)%dinsect_bb_path) // fname

    input_file = jsb_netcdf_open_input(TRIM(fullpath), model%grid_id)
    ptr_2D => input_file%Read_2d(variable_name='dist_fract', fill_array = dinsect_bb_dist_fraction_read(:,:))

    CALL input_file%Close()

    IF (debug_on()) CALL message(TRIM(routine), 'Finishing routine')

  END SUBROUTINE read_bb_dist_fraction

#endif
END MODULE mo_q_dinsect_init
