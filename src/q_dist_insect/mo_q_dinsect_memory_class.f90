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
!>#### declare and define insect disturbance memory class
!>
MODULE mo_q_dinsect_memory_class

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_util,                   ONLY: One_of, int2string
  !USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable
  USE mo_jsb_process_class,      ONLY: Q_DINSECT_
  USE mo_lnd_bgcm_class,         ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
    &                                  LAST_ELEM_ID, get_name_for_element_id, get_short_name_for_element_id, &
    &                                  t_lnd_bgcm
  USE mo_lnd_bgcm_store_class,   ONLY: DINSECT_BGCM_DEADSAPWOOD_ID, DINSECT_BGCM_DEADHEARTWOOD_ID, &
    &                                  DINSECT_BGCM_DEADLEAF_ID, DINSECT_BGCM_VEG_EXTRACT_ID, &
    &                                  DINSECT_BGCM_FAST_LITTER_ID, DINSECT_BGCM_SLOW_LITTER_ID, &
    &                                  DINSECT_BGCM_LARVA_ID, DINSECT_BGCM_DROPPING_ID, &
    &                                  DINSECT_BGCM_SEEDB_FASTL_ID
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_q_dinsect_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 1200 !> TYPE t_jsb_var | 'normal' var: ~ 220 | 5 pools: 270 var | 4 tissues: 24

  ! ======================================================================================================= !
  !> Type definition for dinsect memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_q_dinsect_memory

    !------------------------------------------------------------------------------------------------------ !
    !> QUINCY element-pool structures
    !>
    !TYPE(t_quincy_struct_dinsect_pool)    :: dpool   ! nested pool TYPE t_jsb_n_struct

    !------------------------------------------------------------------------------------------------------ !
    !> bgc_material structures
    !>
    TYPE(t_dinsect_bgcm_main),  POINTER :: dinsect_bgcm_main

    !------------------------------------------------------------------------------------------------------ !
    !------------------------------------------------------------------------------------------------------ !
    !------------------------------------------------------------------------------------------------------ !
    !> Compensatory feeding of defoliator
    !>
    !TYPE(t_jsb_var_real2d)          :: &
    !                                compen_defoli_scale        , &   !< 0~1 0: no compensatory feeding; 1: completely compensatory feeding
    !------------------------------------------------------------------------------------------------------ !
    TYPE(t_jsb_var_real2d)          ::  &
                                    defoliator_effects_triggered, &         !< the flag of calculating timelasting defoliator effects. Currenly, once triggered from 0 to 1, it will keep unchanged.
                                    barkbeetle_effects_triggered, &         !< the flag of calculating timelasting bark beetle effects. Currenly, once triggered from 0 to 1, it will keep unchanged.
                                    defoliated_state, &                    !< the flag for notifying this timestep is under defoliation (1.0_wp) or not (0.0_wp).
                                    is_insect_year, &                      !< the flag for notifying insect disturbance happened on this year (1.0_wp=true) or not (0.0_wp=false).
                                    dinsect_year_since_start, &            !< the number of years since the start of insect disturbance
                                    dinsect_defoliator_resp_c, &           !< the respiration of defoliator [mol C/m2/timestep]
                                    dinsect_defoliator_resp_c13, &         !< the respiration of defoliator C13 [mol C13/m2/timestep]
                                    dinsect_defoliator_resp_c14, &         !< the respiration of defoliator C14 [mol C14/m2/timestep]
                                    dinsect_barkbeetle_labile_c, &         !< the defence used labile C from bark beetle attack [mol C/m2/timestep]
                                    dinsect_barkbeetle_labile_c13, &       !< the defence used labile C13 from bark beetle attack [mol C13/m2/timestep]
                                    dinsect_barkbeetle_labile_c14          !< the defence used labile C14 from bark beetle attack [mol C14/m2/timestep]
    TYPE(t_jsb_var_real2d)          :: &
                                    dinsect_bb_mort_dt, &                  !< Disturbance fraction mortality caused by bark beetle [fraction/timestep]
                                    dinsect_df_mort_dt, &                  !< Disturbance fraction mortality caused by defoliators [fraction/timestep]
                                    dinsect_total_mort_dt                  !< Disturbance fraction mortality caused by all insects disturbances [fraction/timestep]
    TYPE(t_jsb_var_real2d)          :: &
                                    dinsect_pool_total_c                   !< Total carbon from insect pools [mol C/m2/timestep]

  !------------------------------------------------------------------------------------------------------ !

  CONTAINS
    PROCEDURE :: Init => Init_dinsect_memory
  END TYPE t_q_dinsect_memory

  ! ======================================================================================================= !
  !> Type definition for dinsect_bgc_material_elements
  !>
  TYPE, EXTENDS(t_lnd_bgcm) :: t_dinsect_bgcm_with_elements
    TYPE(t_jsb_var_real2d), POINTER :: carbon
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen
    TYPE(t_jsb_var_real2d), POINTER :: phosphorus
    TYPE(t_jsb_var_real2d), POINTER :: carbon13
    TYPE(t_jsb_var_real2d), POINTER :: carbon14
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen15

  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => dinsect_bgcm_with_elements_get_element_name_by_id
    PROCEDURE :: Init                       => dinsect_bgcm_with_elements_init
  END TYPE t_dinsect_bgcm_with_elements

  ! ======================================================================================================= !
  !> Type definition for dinsect_bgc_material_components
  !>
  TYPE, EXTENDS(t_lnd_bgcm) :: t_dinsect_bgcm_with_components
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: leaf
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: fine_root
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: coarse_root
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: sap_wood
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: heart_wood
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: labile
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: reserve
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: fruit
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => dinsect_bgcm_with_components_get_element_name_by_id
    PROCEDURE :: Init                       => dinsect_bgcm_with_components_init
  END TYPE t_dinsect_bgcm_with_components

  ! ======================================================================================================= !
  !> Type definition for dinsect_bgc_material_main
  !>
  !> main structure, added fastly into the "TYPE(t_jsb_memory) :: t_q_dinsect_memory" as "CLASS(t_jsb_pool) :: bgc_material"
  !>
  TYPE, EXTENDS(t_lnd_bgcm) :: t_dinsect_bgcm_main
    !! pools
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: dinsectbdeadsapwood
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: dinsectbdeadheartwood
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: dinsectbdeadleaf
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: dinsectblarva

    !! fluxes
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: dinsectbdropping
    TYPE(t_dinsect_bgcm_with_components),  POINTER :: dinsectbvegextract
    TYPE(t_dinsect_bgcm_with_components),  POINTER :: dinsectbfastlitterfall
    TYPE(t_dinsect_bgcm_with_elements),    POINTER :: dinsectbseedbedfastlitterfall
    TYPE(t_dinsect_bgcm_with_components),  POINTER :: dinsectbslowlitterfall
  CONTAINS
    PROCEDURE :: Get_element_name_by_id     => dinsect_bgcm_main_get_element_name_by_id
    PROCEDURE :: Init                       => dinsect_bgcm_main_init
  END TYPE t_dinsect_bgcm_main

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_memory_class'

CONTAINS
  ! ======================================================================================================= !
  !> Initialize memory (variables) for the process: dist_insect
  !>
  SUBROUTINE Init_dinsect_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    dsl4jsb_Use_config(Q_DINSECT_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_dinsect_memory), INTENT(inout), TARGET  :: mem         !> dinsect memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: lib_id          !< id of primary lct in lctlib
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id

    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),                POINTER :: model                          !< model
    TYPE(t_jsb_grid),                 POINTER :: hgrid                          !< Horizontal grid
    TYPE(t_jsb_vgrid),                POINTER :: surface                        !< Vertical grid
    TYPE(t_jsb_vgrid),                POINTER :: vgrid_soil_w                   !< Vertical grid
    TYPE(t_jsb_vgrid),                POINTER :: vgrid_canopy_q_assimi          !< Vertical grid
    CHARACTER(len=30)                         :: unit_dinsect_pool              !< unit of element var
    CHARACTER(len=30)                         :: unit_dinsect_flux              !< unit of element var
    INTEGER                                   :: elem_idx_map(LAST_ELEM_ID)     !< element mapper ID -> IND
    !TYPE(t_dinsect_bgcm_with_components), POINTER :: dinsect_bgcm_components    !< bgc_material components
    !TYPE(t_dinsect_bgcm_with_elements),   POINTER :: dinsect_bgcm_elements      !< bgc_material elements
    INTEGER                                   :: table                          !< ...
    CHARACTER(len=*), PARAMETER               :: routine = TRIM(modname)//':Init_dinsect_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DINSECT_)
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    surface       => Get_vgrid('surface')
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    vgrid_canopy_q_assimi  => Get_vgrid('q_canopy_layer')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DINSECT_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! init local variables
    unit_dinsect_pool = "mol m-2"
    unit_dinsect_flux = "mol m-2 timestep-1"
    ! ----------------------------------------------------------------------------------------------------- !
    ! add memory only for LAND & PFT tiles
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      ! --------------------------------------------------------------------------------------------------- !
      !>  1.0 bgc_material structures
      !>
      !>
      ! copy elements_index_map from config to local
      elem_idx_map(:) = model%config%elements_index_map(:)
          ! --------------------------------------------------------------------------------------------------- !
      !> create and populate the main bgc_material structure mem%bgc_material that contains all other bgc_material
      !> note: no spaces in 'name' and/or 'shortname'
      !>
      ALLOCATE(t_dinsect_bgcm_main :: mem%bgc_material)
      CALL mem%bgc_material%Init(name = 'dinsect_bgc_material_main', shortname = 'dinsect_bgcm_main', &
        &                        elements_index_map = elem_idx_map(:), l_elements = .FALSE.)
      SELECT TYPE (selected => mem%bgc_material)
      CLASS IS (t_dinsect_bgcm_main)
        mem%dinsect_bgcm_main => selected
      END SELECT


      ! --------------------------------------------------------------------------------------------------- !
      !> Add the info to this memory that this memory contains bgc materials
      !>
      mem%has_bgc_materials = .TRUE.

      ! --------------------------------------------------------------------------------------------------- !
      !> create the bgc_material components structures
      !>
      ! "shortname" is used for model output names
      ! dinsect_bgcm_veg_extract
      mem%dinsect_bgcm_main%dinsectbvegextract => Add_dinsect_bgc_material_with_components(mem%dinsect_bgcm_main, &
        & DINSECT_BGCM_VEG_EXTRACT_ID, name = 'dinsect_bgcm_veg_extract', shortname = 'bgcm_veg_extract',     &
        & unit = unit_dinsect_flux, elements_index_map = elem_idx_map(:))

      ! dinsect_bgcm_fast_litterfall
      mem%dinsect_bgcm_main%dinsectbfastlitterfall => Add_dinsect_bgc_material_with_components(mem%dinsect_bgcm_main, &
        & DINSECT_BGCM_FAST_LITTER_ID, name = 'dinsect_bgcm_fast_litterfall', shortname = 'bgcm_flitterfall',     &
        & unit = unit_dinsect_flux, elements_index_map = elem_idx_map(:))

      ! dinsect_bgcm_slow_litterfall
      mem%dinsect_bgcm_main%dinsectbslowlitterfall => Add_dinsect_bgc_material_with_components(mem%dinsect_bgcm_main, &
        & DINSECT_BGCM_SLOW_LITTER_ID, name = 'dinsect_bgcm_slow_litterfall', shortname = 'bgcm_slitterfall',     &
        & unit = unit_dinsect_flux, elements_index_map = elem_idx_map(:))

      ! --------------------------------------------------------------------------------------------------- !
      !> create the bgc_material elements structures
      !>
      ! "shortname" is used for model output names
      ! dinsect pools and fluxes:
      ! for bark beetles:

      ! dinsect_bgcm_deadsapwood
      mem%dinsect_bgcm_main%dinsectbdeadsapwood => Add_dinsect_bgc_material_with_elements(mem%dinsect_bgcm_main,    &
        & name = 'dinsect_bgcm_deadsapwood', shortname = 'dinsect_bgcm_deadsapwood',               &
        & elements_index_map = elem_idx_map(:), unit = unit_dinsect_pool, bgcm_id = DINSECT_BGCM_DEADSAPWOOD_ID)

      ! dinsect_bgcm_deadheartwood
      mem%dinsect_bgcm_main%dinsectbdeadheartwood => Add_dinsect_bgc_material_with_elements(mem%dinsect_bgcm_main,    &
        & name = 'dinsect_bgcm_deadheartwood', shortname = 'dinsect_bgcm_deadheartwood',               &
        & elements_index_map = elem_idx_map(:), unit = unit_dinsect_pool, bgcm_id = DINSECT_BGCM_DEADHEARTWOOD_ID)

      ! dinsect_bgcm_deadleaf
      mem%dinsect_bgcm_main%dinsectbdeadleaf => Add_dinsect_bgc_material_with_elements(mem%dinsect_bgcm_main,    &
        & name = 'dinsect_bgcm_deadleaf', shortname = 'dinsect_bgcm_deadleaf',               &
        & elements_index_map = elem_idx_map(:), unit = unit_dinsect_pool, bgcm_id = DINSECT_BGCM_DEADLEAF_ID)

      ! for defoliators:
      ! dinsect_bgcm_larva
      mem%dinsect_bgcm_main%dinsectblarva => Add_dinsect_bgc_material_with_elements(mem%dinsect_bgcm_main,    &
        & name = 'dinsect_bgcm_larva', shortname = 'dinsect_bgcm_larva',               &
        & elements_index_map = elem_idx_map(:), unit = unit_dinsect_pool, bgcm_id = DINSECT_BGCM_LARVA_ID)

      ! dinsect_bgcm_dropping
      mem%dinsect_bgcm_main%dinsectbdropping => Add_dinsect_bgc_material_with_elements(mem%dinsect_bgcm_main,    &
        & name = 'dinsect_bgcm_dropping', shortname = 'dinsect_bgcm_dropping',               &
        & elements_index_map = elem_idx_map(:), unit = unit_dinsect_flux, bgcm_id = DINSECT_BGCM_DROPPING_ID)

      ! dinsect_bgcm_seed_bed_fast_litterfall
      mem%dinsect_bgcm_main%dinsectbseedbedfastlitterfall => Add_dinsect_bgc_material_with_elements(mem%dinsect_bgcm_main,    &
        & name = 'dinsect_bgcm_seed_bed_fast_litterfall', shortname = 'bgcm_seed_bed_flitterfall',               &
        & elements_index_map = elem_idx_map(:), unit = unit_dinsect_flux, bgcm_id = DINSECT_BGCM_SEEDB_FASTL_ID)

      ! create list of CHARACTERs with names dinsect_bgcm_main to the "lowest childs" of bgcm structures
      CALL mem%bgc_material%Set_paths()

      ! add variables to bgc_material
      CALL mem%t_jsb_memory_add_bgcm(mem%dinsect_bgcm_main%dinsectbvegextract, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      CALL mem%Add_var(mem%dinsect_bgcm_main%dinsectbfastlitterfall, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      CALL mem%Add_var(mem%dinsect_bgcm_main%dinsectbseedbedfastlitterfall, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      CALL mem%Add_var(mem%dinsect_bgcm_main%dinsectbslowlitterfall, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      CALL mem%Add_var(mem%dinsect_bgcm_main%dinsectbdeadsapwood, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      CALL mem%Add_var(mem%dinsect_bgcm_main%dinsectbdeadheartwood, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      CALL mem%Add_var(mem%dinsect_bgcm_main%dinsectbdeadleaf, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      CALL mem%Add_var(mem%dinsect_bgcm_main%dinsectblarva, &
        & hgrid, surface, vgrid_canopy_q_assimi,                                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      CALL mem%Add_var(mem%dinsect_bgcm_main%dinsectbdropping, &
        & hgrid, surface, vgrid_soil_w,                                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & output_level = BASIC,                                                                        &
        & loutput = .FALSE.,                                                                           &
        & lrestart = .FALSE.,                                                                          &
        & initval_r = 0.0_wp)

      ! --------------------------------------------------------------------------------------------------- !
      !> 2.0 2D & 3D variables
      !>

      CALL mem%Add_var('defoliator_effects_triggered', mem%defoliator_effects_triggered, &
        & hgrid, surface, &
        & t_cf('defoliator_effects_triggered', 'unitless, > 0 means triggered', 'defoliator_effects_triggered'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('barkbeetle_effects_triggered', mem%barkbeetle_effects_triggered, &
        & hgrid, surface, &
        & t_cf('barkbeetle_effects_triggered', 'unitless, > 0 means triggered', 'barkbeetle_effects_triggered'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_bb_mort_dt', mem%dinsect_bb_mort_dt, &
        & hgrid, surface, &
        & t_cf('dinsect_bb_mort_dt', 'fraction ts-1', 'dinsect_bb_mort_dt'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_df_mort_dt', mem%dinsect_df_mort_dt, &
        & hgrid, surface, &
        & t_cf('dinsect_df_mort_dt', 'fraction ts-1', 'dinsect_df_mort_dt'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_total_mort_dt', mem%dinsect_total_mort_dt, &
        & hgrid, surface, &
        & t_cf('dinsect_total_mort_dt', 'fraction ts-1', 'dinsect_total_mort_dt'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('defoliated_state', mem%defoliated_state, &
        & hgrid, surface, &
        & t_cf('defoliated_state', 'unitless, > 0 means defoliated', 'defoliated_state'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('is_insect_year', mem%is_insect_year, &
        & hgrid, surface, &
        & t_cf('is_insect_year', 'unitless, > 0 means is insect year', 'is_insect_year'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_year_since_start', mem%dinsect_year_since_start, &
        & hgrid, surface, &
        & t_cf('dinsect_year_since_start', 'years', 'dinsect_year_since_start'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_defoliator_resp_c', mem%dinsect_defoliator_resp_c, &
        & hgrid, surface, &
        & t_cf('dinsect_defoliator_resp_c', 'mol m-2 timestep-1', 'dinsect_defoliator_resp_c'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_defoliator_resp_c13', mem%dinsect_defoliator_resp_c13, &
        & hgrid, surface, &
        & t_cf('dinsect_defoliator_resp_c13', 'mol m-2 timestep-1', 'dinsect_defoliator_resp_c13'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_defoliator_resp_c14', mem%dinsect_defoliator_resp_c14, &
        & hgrid, surface, &
        & t_cf('dinsect_defoliator_resp_c14', 'mol m-2 timestep-1', 'dinsect_defoliator_resp_c14'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_barkbeetle_labile_c', mem%dinsect_barkbeetle_labile_c, &
        & hgrid, surface, &
        & t_cf('dinsect_barkbeetle_labile_c', 'mol m-2 timestep-1', 'dinsect_barkbeetle_labile_c'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_barkbeetle_labile_c13', mem%dinsect_barkbeetle_labile_c13, &
        & hgrid, surface, &
        & t_cf('dinsect_barkbeetle_labile_c13', 'mol m-2 timestep-1', 'dinsect_barkbeetle_labile_c13'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_barkbeetle_labile_c14', mem%dinsect_barkbeetle_labile_c14, &
        & hgrid, surface, &
        & t_cf('dinsect_barkbeetle_labile_c14', 'mol m-2 timestep-1', 'dinsect_barkbeetle_labile_c14'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('dinsect_pool_total_c', mem%dinsect_pool_total_c, &
        & hgrid, surface, &
        & t_cf('dinsect_pool_total_c', 'mol C m-2 timestep-1', 'Total Carbon from insect pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .FALSE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

  ENDIF
END SUBROUTINE Init_dinsect_memory

! ======================================================================================================= !
!>
!> return name of element variable of bgc_material by its ENUM ID
!> helper function for t_dinsect_bgcm_with_elements
!> "bgcm_object%Get_element_name_by_id" => dinsect_bgcm_with_elements_get_element_name_by_id
!>
FUNCTION dinsect_bgcm_with_elements_get_element_name_by_id(this, id) RESULT(name)
  CLASS(t_dinsect_bgcm_with_elements), INTENT(in)  :: this
  INTEGER,                         INTENT(in)  :: id
  CHARACTER(LEN=:),                ALLOCATABLE :: name

  ! return name of element variable only if this bgc_material object contains elements
  IF (this%contains_elements) THEN
    name = get_name_for_element_id(id)
  ELSE
    name='no_elements'
  END IF
END FUNCTION dinsect_bgcm_with_elements_get_element_name_by_id

! ======================================================================================================= !
!>
!> return name of element variable of bgc_material by its ENUM ID
!> helper function for t_dinsect_bgcm_with_components
!> "bgcm_object%Get_element_name_by_id" => dinsect_bgcm_with_components_get_element_name_by_id
!>
FUNCTION dinsect_bgcm_with_components_get_element_name_by_id(this, id) RESULT(name)
  CLASS(t_dinsect_bgcm_with_components), INTENT(in)  :: this
  INTEGER,                           INTENT(in)  :: id
  CHARACTER(LEN=:),                  ALLOCATABLE :: name

  ! return name of element variable only if this bgc_material object contains elements
  IF (this%contains_elements) THEN
    name = get_name_for_element_id(id)
  ELSE
    name='no_elements'
  END IF
END FUNCTION dinsect_bgcm_with_components_get_element_name_by_id

! ======================================================================================================= !
!>
!> return name of element variable of bgc_material by its ENUM ID
!> helper function for t_dinsect_bgcm_with_components
!> "bgcm_object%Get_element_name_by_id" => dinsect_bgcm_main_get_element_name_by_id
!>
FUNCTION dinsect_bgcm_main_get_element_name_by_id(this, id) RESULT(name)
  CLASS(t_dinsect_bgcm_main), INTENT(in)  :: this
  INTEGER,                INTENT(in)  :: id
  CHARACTER(LEN=:),       ALLOCATABLE :: name

  ! return name of element variable only if this bgc_material object contains elements
  IF (this%contains_elements) THEN
    name = get_name_for_element_id(id)
  ELSE
    name='no_elements'
  END IF
END FUNCTION dinsect_bgcm_main_get_element_name_by_id

! ======================================================================================================= !
!>
!> initializing the t_dinsect_bgcm_with_elements object
!> "bgcm_object%Init()"
!> by default it does contain element variables !
!>
SUBROUTINE dinsect_bgcm_with_elements_init(this, name, shortname, elements_index_map, &
  &                                    l_elements, element_list, element_unit)
  ! ----------------------------------------------------------------------------------------------------- !
  CLASS(t_dinsect_bgcm_with_elements),    INTENT(inout) :: this               ! object of TYPE t_dinsect_elements_pool
  CHARACTER(LEN=*),                   INTENT(in)    :: name                   ! pool name       (No spaces !)
  CHARACTER(LEN=*),                   INTENT(in)    :: shortname              ! pool short name (No spaces !)
  INTEGER,                            INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
  LOGICAL,          OPTIONAL,         INTENT(in)    :: l_elements             ! IF pool has elements
  INTEGER,          OPTIONAL,         INTENT(in)    :: element_list(:) ! @TODO remove | list of elements; used if only selected elements should be added
  CHARACTER(LEN=*), OPTIONAL,         INTENT(in)    :: element_unit           ! unit of elements
  ! ----------------------------------------------------------------------------------------------------- !
  INTEGER                     :: element_id
  CHARACTER(len=*), PARAMETER :: routine = modname//':dinsect_bgcm_with_elements_init'
  ! ----------------------------------------------------------------------------------------------------- !

  this%name              = name
  this%shortname         = shortname
  this%contains_elements = l_elements
  this%contains_bgcm     = .FALSE.

  IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')
  this%element_unit = element_unit

  ! add element variables
  ! add element variables
  DO element_id = 1, SIZE(elements_index_map)
    IF (elements_index_map(element_id) > 0) THEN
      CALL this%Add_element(get_name_for_element_id(element_id), get_short_name_for_element_id(element_id), element_id, dim=2)

      SELECT TYPE (selector => this%element_list(elements_index_map(element_id))%p)
      CLASS IS (t_jsb_var_real2d)
        SELECT CASE(element_id)
        CASE(ELEM_C_ID)
          this%carbon => selector
        CASE(ELEM_N_ID)
          this%nitrogen => selector
        CASE(ELEM_P_ID)
          this%phosphorus => selector
        CASE(ELEM_C13_ID)
          this%carbon13 => selector
        CASE(ELEM_C14_ID)
          this%carbon14 => selector
        CASE(ELEM_N15_ID)
          this%nitrogen15 => selector
        CASE DEFAULT
          CALL finish(routine, 'Invalid element id'//int2string(element_id))
        END SELECT
      END SELECT
    ENDIF
  ENDDO

END SUBROUTINE dinsect_bgcm_with_elements_init

! ======================================================================================================= !
!>
!> initializing the t_dinsect_bgcm_with_components object
!> "bgcm_object%Init()"
!> by default it does not contain element variables (only t_dinsect_bgcm_with_elements do)
!>
SUBROUTINE dinsect_bgcm_with_components_init(this, name, shortname, elements_index_map, &
  &                                      l_elements, element_list, element_unit)
  CLASS(t_dinsect_bgcm_with_components),    INTENT(inout) :: this               ! object of TYPE t_dinsect_components_pool
  CHARACTER(LEN=*),                     INTENT(in)    :: name                   ! pool name       (No spaces !)
  CHARACTER(LEN=*),                     INTENT(in)    :: shortname              ! pool short name (No spaces !)
  INTEGER,                              INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
  LOGICAL,          OPTIONAL,           INTENT(in)    :: l_elements             ! IF pool has elements
  INTEGER,          OPTIONAL,           INTENT(in)    :: element_list(:)        ! list of elements; used if only selected elements should be added
  CHARACTER(LEN=*), OPTIONAL,           INTENT(in)    :: element_unit           ! unit of elements
  ! ----------------------------------------------------------------------------------------------------- !
  CHARACTER(len=*), PARAMETER :: routine = modname//':dinsect_bgcm_with_components_init'

  this%name              = name
  this%shortname         = shortname
  this%contains_elements = l_elements
  this%contains_bgcm     = .TRUE.

  IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')
  this%element_unit = element_unit
END SUBROUTINE dinsect_bgcm_with_components_init

! ======================================================================================================= !
!>
!> initializing the t_dinsect_bgcm_main object
!> "bgcm_object%Init()"
!> by default it does not contain element variables (only t_dinsect_bgcm_with_elements do)
!>
SUBROUTINE dinsect_bgcm_main_init(this, name, shortname, elements_index_map, &
  &                           l_elements, element_list, element_unit)
  ! ----------------------------------------------------------------------------------------------------- !
  CLASS(t_dinsect_bgcm_main),     INTENT(inout) :: this               ! object of TYPE t_dinsect_components_pool
  CHARACTER(LEN=*),               INTENT(in)    :: name                   ! pool name       (No spaces !)
  CHARACTER(LEN=*),               INTENT(in)    :: shortname              ! pool short name (No spaces !)
  INTEGER,                        INTENT(in)    :: elements_index_map(:)  ! elements ID -> IND
  LOGICAL,          OPTIONAL,     INTENT(in)    :: l_elements             ! IF pool has elements
  INTEGER,          OPTIONAL,     INTENT(in)    :: element_list(:)        ! list of elements; used if only selected elements should be added
  CHARACTER(LEN=*), OPTIONAL,     INTENT(in)    :: element_unit           ! unit of elements
  ! ----------------------------------------------------------------------------------------------------- !
  CHARACTER(len=*), PARAMETER :: routine = modname//':dinsect_bgcm_main_init'

  this%name              = name
  this%shortname         = shortname
  this%contains_elements = l_elements
  this%contains_bgcm     = .TRUE.
END SUBROUTINE dinsect_bgcm_main_init

! ======================================================================================================= !
! >
! > create and init an object of the Q_DINSECT_ components bgcm "t_insect_bgcm_with_components"
! > used for DINSECT_ memory init
! >   init and add all its objects of the "t_jsb_pool" with "pool_object%Init()"
! >   "object => Create_dinsect_bgc_material_components()"
! >
FUNCTION Add_dinsect_bgc_material_with_components(pool, bgcm_id, name, shortname, unit, elements_index_map) &
    & RESULT(dinsect_bgcm_with_components)
  ! ----------------------------------------------------------------------------------------------------- !
  CLASS(t_lnd_bgcm), TARGET,  INTENT(INOUT)  :: pool                  !< pool to which this pool with component should be added
  INTEGER,                    INTENT(in)     :: bgcm_id               !< id of this dinsect bgcm with components
  CHARACTER(LEN=*),           INTENT(in)     :: name                  !< name of this dinsect bgcm with components
  CHARACTER(LEN=*),           INTENT(in)     :: shortname             !< shortname of this dinsect bgcm with components
  CHARACTER(LEN=*),           INTENT(in)     :: unit                  !< unit of this dinsect bgcm with components
  INTEGER,                    INTENT(in)     :: elements_index_map(:) !< elements ID -> IND
  ! ----------------------------------------------------------------------------------------------------- !
  TYPE(t_dinsect_bgcm_with_components), POINTER :: dinsect_bgcm_with_components !< the new dinsect bgcm with components
  ! ----------------------------------------------------------------------------------------------------- !

  ! bgcm_components init
  ALLOCATE(dinsect_bgcm_with_components)
  CALL dinsect_bgcm_with_components%Init(name, shortname, elements_index_map(:), l_elements=.FALSE., element_unit=unit)

  ! Add the components (leaf, fine_root, etc) using Add_dinsect_bgc_material_with_elements
  ! leaf
  dinsect_bgcm_with_components%leaf => Add_dinsect_bgc_material_with_elements(dinsect_bgcm_with_components, &
    & name = 'leaf_bgcm', shortname = 'leaf', elements_index_map = elements_index_map(:), unit = unit)
  ! fine_root
  dinsect_bgcm_with_components%fine_root => Add_dinsect_bgc_material_with_elements(dinsect_bgcm_with_components, &
    & name = 'fine_root_bgcm', shortname = 'fine_root', elements_index_map = elements_index_map(:), unit = unit)
  ! coarse_root
  dinsect_bgcm_with_components%coarse_root => Add_dinsect_bgc_material_with_elements(dinsect_bgcm_with_components, &
    & name = 'coarse_root_bgcm', shortname = 'coarse_root', elements_index_map = elements_index_map(:), unit = unit)
  ! sap_wood
  dinsect_bgcm_with_components%sap_wood => Add_dinsect_bgc_material_with_elements(dinsect_bgcm_with_components, &
    & name = 'sap_wood_bgcm', shortname = 'sap_wood', elements_index_map = elements_index_map(:), unit = unit)
  ! heart_wood
  dinsect_bgcm_with_components%heart_wood => Add_dinsect_bgc_material_with_elements(dinsect_bgcm_with_components, &
    & name = 'heart_wood_bgcm', shortname = 'heart_wood', elements_index_map = elements_index_map(:), unit = unit)
  ! labile
  dinsect_bgcm_with_components%labile => Add_dinsect_bgc_material_with_elements(dinsect_bgcm_with_components, &
    & name = 'labile_bgcm', shortname = 'labile', elements_index_map = elements_index_map(:), unit = unit)
  ! reserve
  dinsect_bgcm_with_components%reserve => Add_dinsect_bgc_material_with_elements(dinsect_bgcm_with_components, &
    & name = 'reserve_bgcm', shortname = 'reserve', elements_index_map = elements_index_map(:), unit = unit)
  ! fruit
  dinsect_bgcm_with_components%fruit => Add_dinsect_bgc_material_with_elements(dinsect_bgcm_with_components, &
    & name = 'fruit_bgcm', shortname = 'fruit', elements_index_map = elements_index_map(:), unit = unit)

  ! Add the bookkeeping for this dinsect bgcm with components
  CALL pool%Add_bgcm(dinsect_bgcm_with_components, bgcm_id)

END FUNCTION Add_dinsect_bgc_material_with_components

FUNCTION Add_dinsect_bgc_material_with_elements(pool, name, shortname, elements_index_map, unit, bgcm_id) &
    & RESULT(dinsect_bgcm_with_elements)
  ! ----------------------------------------------------------------------------------------------------- !
  CLASS(t_lnd_bgcm), TARGET,  INTENT(INOUT) :: pool                   !< pool to which this pool with elements should be added
  CHARACTER(LEN=*),           INTENT(in)    :: name                   !< name of this dinsect bgcm with elements
  CHARACTER(LEN=*),           INTENT(in)    :: shortname              !< shortname of this dinsect bgcm with elements
  INTEGER,                    INTENT(in)    :: elements_index_map(:)  !< mapping elements ID -> IND
  CHARACTER(LEN=*),           INTENT(in)    :: unit                   !< unit of this dinsect bgcm with elements
  INTEGER,          OPTIONAL, INTENT(in)    :: bgcm_id                !< optional: id of this bgcm with elements
  ! ----------------------------------------------------------------------------------------------------- !
  TYPE(t_dinsect_bgcm_with_elements), POINTER :: dinsect_bgcm_with_elements !< the new dinsect bgcm with elements
  ! ----------------------------------------------------------------------------------------------------- !

  ALLOCATE(dinsect_bgcm_with_elements)
  CALL dinsect_bgcm_with_elements%Init(name = name, shortname = shortname, &
    & elements_index_map = elements_index_map(:), l_elements = .TRUE., element_unit = unit)
  CALL pool%Add_bgcm(dinsect_bgcm_with_elements, bgcm_id)

END FUNCTION Add_dinsect_bgc_material_with_elements

END MODULE mo_q_dinsect_memory_class
