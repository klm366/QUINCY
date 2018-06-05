!>
!! @brief Contains methods for JSBACH initializations.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                          by Reiner Schnur (2013-04-04)
!! Some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_model_init
#ifndef __NO_JSBACH__

  USE mo_jsb_base,         ONLY: jsbach
  USE mo_jsb_control,      ONLY: get_no_of_models, debug_on
  USE mo_jsb_model_class,  ONLY: t_jsb_model

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message, finish
  USE mo_jsb_parallel,     ONLY: Get_omp_thread, Is_omp_inside_serial

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_model

  INTERFACE init_model
    MODULE PROCEDURE init_model_standalone
    MODULE PROCEDURE init_model_patch
    MODULE PROCEDURE init_model_decomposition
  END INTERFACE init_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_init'

CONTAINS

  SUBROUTINE init_model_standalone(nproca, nprocb, npedim, model_id, model_name)
#ifndef __ICON__
    !>>> TODO: echam
    USE mo_gaussgrid,     ONLY : inigau
    USE mo_echam_convect_tables, ONLY : init_convect_tables
    USE mo_netCDF,      ONLY: add_dim
    !<<<
#endif
    USE mo_jsb_domain_iface, ONLY: get_nlon_nlat_from_file, get_lons_and_lats_from_file, jsbalone_init_decomp, pe_decomposed
    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,         ONLY: Register_grid
    USE mo_jsb_model_class,  ONLY: print_model


    INTEGER,                    INTENT(in) :: model_id, nproca, nprocb, npedim
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: model_name

    TYPE(t_jsb_model),  POINTER  :: model
    TYPE(t_jsb_grid),   POINTER  :: grid
    TYPE(pe_decomposed) :: decomposition

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_standalone'

    INTEGER :: nlon, nlat

    model => jsbach%models(model_id)%m

    IF (PRESENT(model_name)) THEN
      IF (ADJUSTL(TRIM(model_name)) /= ADJUSTL(TRIM(model%name))) &
        CALL finish(TRIM(routine), 'Model ID and name are inconsistent with namelist config.')
    END IF

    CALL get_nlon_nlat_from_file(model%config%grid_filename, nlon, nlat)
    CALL jsbalone_init_decomp(nlon, nlat, nproca, nprocb, npedim, decomposition)

#ifndef __ICON__
    CALL inigau               ! beside others required for get_area, i.e. required before new_grid
    CALL init_convect_tables

    CALL add_dim("lon", nlon, "longitude", "degrees_east" ) !required to write restart files
    CALL add_dim("lat", nlat, "latitude", "degrees_north")

#endif

    !needs to be done after the decomposition is initialised; required for get_lon and get_lat, especially for site-level sims
#ifndef __ICON__
    CALL get_lons_and_lats_from_file(model%config%grid_filename, decomposition%nproma, decomposition%ngpblks)
#else
    CALL finish(TRIM(routine),'Standalone model currently only works with ECHAM')
#endif

    grid => new_grid(decomp=decomposition)

    CALL register_grid(grid)
    model%grid_id = grid%id

    CALL init_model_common(model)

    NULLIFY(model)

  END SUBROUTINE init_model_standalone

  SUBROUTINE init_model_patch( model_id, patch)

    USE mo_jsb_domain_iface, ONLY: t_patch
    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,         ONLY: Register_grid

    INTEGER,                    INTENT(in) :: model_id
    TYPE(t_patch),              INTENT(in) :: patch

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_patch'

    IF (patch%id < 0) &
      CALL finish(TRIM(routine), 'Improper call with invalid patch')

    IF (model_id > get_no_of_models()) &
      CALL finish(TRIM(routine), 'Host request to more jsbach models than are defined')

    model => jsbach%models(model_id)%m

    grid => new_grid(patch=patch)
    CALL Register_grid(grid)
    model%grid_id = grid%id

    CALL init_model_common(model)

    NULLIFY(model)

  END SUBROUTINE init_model_patch

  SUBROUTINE init_model_decomposition( model_id, decomposition)

    USE mo_jsb_domain_iface, ONLY: pe_decomposed
    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,         ONLY: Register_grid

    INTEGER,               INTENT(in) :: model_id
    TYPE(pe_decomposed),   INTENT(in) :: decomposition

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_decomposition'

    IF (decomposition%nlon < 0) &
      CALL finish(TRIM(routine), 'Improper call with invalid decomposition')

    IF (model_id > get_no_of_models()) &
      CALL finish(TRIM(routine), 'Host request to more jsbach models than are defined')

    model => jsbach%models(model_id)%m

    grid => new_grid(decomp=decomposition)
    CALL Register_grid(grid)
    model%grid_id = grid%id

    CALL init_model_common(model)

    NULLIFY(model)

  END SUBROUTINE init_model_decomposition

  SUBROUTINE init_model_common(model)

    USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_, AGGREGATE_
    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process
    USE mo_jsb_model_class,     ONLY: print_model
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_lct_class,       ONLY: Contains_lct
    USE mo_jsb_lctlib_class,    ONLY: Read_lctlib
    USE mo_jsb_model_usecases,  ONLY: init_usecase
    USE mo_seb_init,            ONLY: seb_init
    USE mo_turb_init,           ONLY: turb_init
    USE mo_sse_init,            ONLY: sse_init
    USE mo_hydro_init,          ONLY: hydro_init
    USE mo_rad_init,            ONLY: rad_init
    !USE mo_assimi_init,         ONLY: assimi_init
    USE mo_pheno_init,          ONLY: pheno_init
    USE mo_carbon_init,         ONLY: carbon_init
    USE mo_fage_init,           ONLY: fage_init
    USE mo_fm_init,             ONLY: fm_init
    USE mo_harv_init,           ONLY: harv_init
    USE mo_fuel_init,           ONLY: fuel_init
    USE mo_disturb_init,        ONLY: disturb_init
#ifndef __NO_JSBACH_HD__
    USE mo_hd_init,             ONLY: hd_init_ic, hd_init_bc
#endif
    USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,            ONLY: Register_vgrid
    USE mo_jsb_io,              ONLY: ZAXIS_SURFACE, cdiDefMissval, missval
    USE mo_util,                ONLY: one_of

    dsl4jsb_Use_processes SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_, CARBON_, DISTURB_, FUEL_, FAGE_, FM_, HARV_
    dsl4jsb_Use_config(PHENO_)
    dsl4jsb_Use_config(SSE_)
    dsl4jsb_Use_config(HYDRO_)
    dsl4jsb_Use_config(RAD_)
    dsl4jsb_Use_config(CARBON_)
    dsl4jsb_Use_config(DISTURB_)
    dsl4jsb_Use_config(FUEL_)
    dsl4jsb_Use_config(FAGE_)
    dsl4jsb_Use_config(FM_)
    dsl4jsb_Use_config(HARV_)

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_config(DISTURB_)
    dsl4jsb_Def_config(FUEL_)
    dsl4jsb_Def_config(FM_)
    dsl4jsb_Def_config(HARV_)
    dsl4jsb_Def_config(FAGE_)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile_abstract), POINTER :: tile, child
    TYPE(t_jsb_vgrid),          POINTER :: surface                      ! Vertical grid
    INTEGER                             :: iproc, ilct, jlct, i, no_children
    REAL(wp), ALLOCATABLE               :: fract_tmp(:,:)
    INTEGER                             :: no_omp_thread

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_common'

    no_omp_thread = Get_omp_thread()

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    CALL model%Set_mode(model%config%hsm_mode)

    CALL cdiDefMissval(missval)

    ALLOCATE(model%processes(max_no_of_processes))
    DO iproc=1,max_no_of_processes
      model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename))
    END DO

    CALL model%Configure_processes()

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(RAD_)
    dsl4jsb_Get_config(CARBON_)
    dsl4jsb_Get_config(FUEL_)
    dsl4jsb_Get_config(DISTURB_)
    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_config(FM_)
    dsl4jsb_Get_config(HARV_)

    ! Crosscheck process configs
    IF (model%config%l_compat401) THEN
      CALL message(TRIM(routine), 'Using settings compatible to jsbach4.01:')
      CALL message('    ', '... using soil heat capacity from mineral soil map')
      CALL message('    ', '... using soil heat conductivity from mineral soil map')
      CALL message('    ', '... using fixed snow heat capacity and snow heat conductivity')
      CALL message('    ', '... not using multi-layer snow model')
      CALL message('    ', '... not using dynamical snow density and thickness')
      CALL message('    ', '... not using frozen ground')
      CALL message('    ', '... not using organic soil fractions')
      dsl4jsb_Config(SSE_)%l_heat_cap_dyn = .FALSE.
      dsl4jsb_Config(SSE_)%l_heat_cond_dyn = .FALSE.
      dsl4jsb_Config(SSE_)%l_heat_cap_map = .FALSE.
      dsl4jsb_Config(SSE_)%l_heat_cond_map = .FALSE.
      dsl4jsb_Config(SSE_)%l_snow = .FALSE.
      dsl4jsb_Config(SSE_)%l_dynsnow = .FALSE.
      dsl4jsb_Config(SSE_)%l_freeze = .FALSE.
      dsl4jsb_Config(HYDRO_)%l_organic = .FALSE.
    END IF

    IF (dsl4jsb_Config(DISTURB_)%active .AND. .NOT. dsl4jsb_Config(FUEL_)%active) THEN
      dsl4jsb_Config(FUEL_)%active = .TRUE.
      CALL message(TRIM(routine), 'Switching on fuel process since disturbance process is active')
    END IF
    IF (dsl4jsb_Config(DISTURB_)%active .AND. .NOT. dsl4jsb_Config(CARBON_)%active) THEN
      dsl4jsb_Config(CARBON_)%active = .TRUE.
      CALL message(TRIM(routine), 'Switching on carbon process since disturbance process is active')
    END IF

    surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
    CALL register_vgrid(surface)

    model%lctlib => Read_lctlib('lctlib_nlct21.def')
    CALL message(TRIM(routine), 'Imported lctlib from '//'lctlib_nlct21.def')

    IF (model%config%usecase /= '') THEN
      ! Shortcut model configuration through usecase
      CALL init_usecase(model)
    ELSE
      ! Explicit model configuration through namelist
      ! TBD
    END IF

    !Aging happens once at the beginning of each year at the start of a timestep. 
    IF (model%processes(FAGE_)%p%config%active) THEN
       CALL model%Queue_task(FAGE_, "calculate_ageing")
       CALL model%Queue_task(CARBON_, "apply_age_cohort_changes")
    END IF

    !Forest Management happens once at the beginning of each year at the start of a timestep - after ageing 
    !NOTE: harvest_age_cohorts ASSUMES that it runs after calculate_ageing, and thus that the youngest age cohort is empty! 
    IF (model%processes(FM_)%p%config%active) THEN
      CALL model%Queue_task(FM_,     "forest_management")
      CALL model%Queue_task(HARV_,   "forest_harvest")
      IF (model%processes(FAGE_)%p%config%active) THEN
        CALL model%Queue_task(FAGE_,   "harvest_age_cohorts")
        CALL model%Queue_task(CARBON_, "apply_age_cohort_changes")
      ENDIF
    ENDIF

    ! Temporary for testing ... should be in separate subroutine
    CALL model%Queue_task(RAD_,    "surface_radiation")
    CALL model%Queue_task(PHENO_,  "phenology")           ! JN moved phenology up - for allom it needs to happen
                                                          ! before C_NPP_pot_allocation TODO-JN
                      ! before it was between humidity_scaling and foliage_projected_cover
                      ! + prev comment: can be placed between the task "surface_radiation" and "asselin_filter"
    CALL model%Queue_task(HYDRO_,  "snow_and_skin_fraction")
    CALL model%Queue_task(HYDRO_,  "soil_properties")
    CALL model%Queue_task(SEB_,    "surface_energy")
    CALL model%Queue_task(HYDRO_,  "evaporation")
    CALL model%Queue_task(SEB_,    "surface_fluxes")
    IF (model%processes(ASSIMI_)%p%config%active) THEN
      CALL model%Queue_task(RAD_,    "radiation_par")
    END IF
    CALL model%Queue_task(HYDRO_,  "surface_hydrology")
    CALL model%Queue_task(HYDRO_,  "snow_and_skin_fraction")
    CALL model%Queue_task(SSE_,    "soil_properties")
    CALL model%Queue_task(SSE_,    "soil_temperature")
    CALL model%Queue_task(HYDRO_,  "soil_hydrology")
#ifndef __NO_JSBACH_HD__
    CALL model%Queue_task(HD_,     "hd_local_outflow")
#endif
    CALL model%Queue_task(SEB_,    "snowmelt_correction")
    IF (model%processes(ASSIMI_)%p%config%active) THEN
      ! Normally called before
      ! assimilation -> update_assimilation: -> calc_assimilation_waterlimited
      ! and before
      ! canopy_cond_unstressed -> update_canopy_conductance_unstressed: 
      !                                    -> calc_assimilation_waterunlimited
      CALL model%Queue_task(ASSIMI_, "assimilation_scaling_factors")
    END IF
    IF (model%processes(ASSIMI_)%p%config%active) THEN
      CALL model%Queue_task(ASSIMI_, "canopy_cond_unstressed") ! -> calc_assimilation_waterunlimited
    ELSE
      CALL model%Queue_task(HYDRO_,  "canopy_cond_unstressed")
    END IF
    CALL model%Queue_task(HYDRO_,  "canopy_cond_stressed") ! Should be called with and without the ASSIMI process,
                                                           ! and before assimilation
    ! CALL model%Queue_task(HYDRO_,  "snow_and_skin_fraction")
    IF (model%processes(ASSIMI_)%p%config%active) THEN
      CALL model%Queue_task(ASSIMI_, "assimilation") ! -> calc_assimilation_waterlimited
                                                     ! should be called after the soil hydrology is calculated
    END IF

    IF (model%processes(FUEL_)%p%config%active) THEN
      CALL model%Queue_task(FUEL_, "fuel")
    END IF
    IF (model%processes(DISTURB_)%p%config%active) THEN
      CALL model%Queue_task(DISTURB_, "natural_disturbances")

      ! in case that forest age is switched on, i.e. working with age cohorts,
      ! disturbances need to be handled on the parent of the age cohorts!
      IF (model%processes(FAGE_)%p%config%active) THEN
         CALL model%Queue_task(FAGE_, "nat_dist_age_cohorts")
         CALL model%Queue_task(CARBON_, "apply_age_cohort_changes")
      END IF
    END IF

    ! JN: changed order npp and disturbances because disturbances might change area 
    !     and by swapping dist and carbon less variables have to be aggregated upon area shift
    IF (model%processes(CARBON_)%p%config%active) THEN
      CALL model%Queue_task(CARBON_, "C_NPP_pot_allocation")
    END IF

    CALL model%Queue_task(TURB_,   "humidity_scaling")

    CALL model%Queue_task(PHENO_,  "foliage_projected_cover")
    CALL model%Queue_task(TURB_,   "surface_roughness")
    ! Note: albedo is updated at the end of the time step, one time step before the radiation is updated in the atmosphere
    CALL model%Queue_task(RAD_,    "albedo")
    CALL model%Queue_task(HYDRO_,  "water_balance")
    CALL model%Queue_task(SEB_,    "asselin_filter")      ! Should always be the last task
    IF (debug_on()) THEN
      CALL print_model(model)
      CALL model%Print_actions()
    END IF

    ! At this point, the model configuration and tile structure is set up.

    ! We continue with initializing the memory structure and different model processes.

    ! Iterate over tiles and call initialization if necessary
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        tile => model%Goto_next_tile()
        CYCLE
      END IF

!TODO-JN: l_fixed_fractions = .FALSE.!
      ! Initialize tile memory and read tile fractions
      ! @todo: once we have land use change or dynveg, l_fixed_fractions needs to be adapted
      CALL tile%Init(TRIM(model%shortname), '', TRIM(tile%name), l_fixed_fractions=.FALSE., &
        & in_var_groups=one_of(TRIM(tile%name), model%config%output_tiles) > 0)

      ! TBD: Loop over processes instead of explicit
      IF (one_of(tile%process_action(RAD_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL rad_init(tile)
      END IF

      IF (one_of(tile%process_action(SEB_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL seb_init(tile)
      END IF

      IF (one_of(tile%process_action(TURB_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL turb_init(tile)
      END IF

      ! Initialize PHENO_ before HYDRO_
      IF (one_of(tile%process_action(PHENO_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL pheno_init(tile)
      END IF

      ! HYDRO_ needs to be initialized before SSE_
      IF (one_of(tile%process_action(HYDRO_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL hydro_init(tile)
      END IF

      IF (one_of(tile%process_action(SSE_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL sse_init(tile)
      END IF

      ! IF (one_of(tile%process_action(ASSIMI_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
      !   CALL assimi_init_bc(tile)
      !   CALL assimi_init_ic(tile)
      ! END IF

      IF (one_of(tile%process_action(CARBON_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL carbon_init(tile)
      END IF

      IF (one_of(tile%process_action(FUEL_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL fuel_init(tile)
      END IF

      ! R: in the moment we do not need any initialization for disturbance:
      !IF (one_of(tile%process_action(DISTURB_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
      !  CALL disturb_init(tile)
      !END IF

      IF (one_of(tile%process_action(FAGE_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL fage_init(tile)
      END IF

      IF (one_of(tile%process_action(FM_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL fm_init(tile)
      END IF

      IF (one_of(tile%process_action(HARV_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL harv_init(tile)
      END IF

#ifndef __NO_JSBACH_HD__
      IF (one_of(tile%process_action(HD_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL hd_init_bc(tile)
        CALL hd_init_ic(tile)
      END IF
#endif

      tile => model%Goto_next_tile()

    END DO

    ! Iterate over tiles again and set lct fractions (after ALL tile fractions have been set above)
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    ALLOCATE(fract_tmp(SIZE(tile%lcts(1)%fract,1), SIZE(tile%lcts(1)%fract,2)))

    DO WHILE (ASSOCIATED(tile))

      IF (.NOT. tile%visited(no_omp_thread) .AND. tile%Has_children()) THEN
        tile => model%Goto_next_tile()
        CYCLE
      END IF

      ! tile is a leaf or it has children and we're on the way up
      DO ilct=1,SIZE(tile%lcts)
        IF (tile%Has_children()) THEN
          tile%lcts(ilct)%fract(:,:) = 0._wp
          no_children = tile%Get_no_of_children()
          child => tile%Get_first_child_tile()
          DO i=1,no_children
            IF (Contains_lct(child%lcts, tile%lcts(ilct)%id)) THEN
              DO jlct=1,SIZE(child%lcts)
                IF (child%lcts(jlct)%id == tile%lcts(ilct)%id) THEN
                  CALL child%Get_fraction(fract=fract_tmp(:,:))
                  IF (jlct > 1) THEN
                    fract_tmp(:,:) = fract_tmp(:,:) * child%lcts(jlct)%fract(:,:)
                  END IF
                  tile%lcts(ilct)%fract(:,:) = tile%lcts(ilct)%fract(:,:) + fract_tmp(:,:)
                  EXIT
                END IF
              END DO
            END IF
            child => child%Get_next_sibling_tile()
          END DO
        ELSE
          IF (ilct > 1) CALL finish(TRIM(routine), 'Leave tiles should only have one LCT')
          tile%lcts(1)%fract(:,:) = 1._wp        ! Leaves have only one LCT
        END IF
      END DO

      tile => model%Goto_next_tile()

    END DO

    DEALLOCATE(fract_tmp)

    ! Iterate over tiles again and initialize hierarchical structure for memories and variables
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (.NOT. tile%visited(no_omp_thread) .AND. tile%Has_children()) THEN
        tile => model%Goto_next_tile()
        CYCLE
      END IF

      CALL tile%Init_vars()

      tile => model%Goto_next_tile()

    END DO

    ! Iterate over tiles and crossckeck tile structure
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        tile => model%Goto_next_tile()
        CYCLE
      END IF

      IF (one_of(tile%process_action(PHENO_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        dsl4jsb_Get_config(PHENO_)
        IF (      tile%is_vegetation &
          & .AND. TRIM(dsl4jsb_Config(PHENO_)%scheme) /= 'climatology' &
          & .AND. tile%lcts(1)%lib_id == 0) &
            & CALL finish(TRIM(routine), 'Must use phenology climatology on general, non-PFT VEG tile '//TRIM(tile%name))
      END IF

      tile => model%Goto_next_tile()

    END DO
    !
    !
    IF (debug_on()) THEN
      ! Iterate over tiles again and print tile information (after ALL tiles have been initialized above)
      CALL model%Reset_tiles()
      tile => model%Get_top_tile()
      IF (.NOT. ASSOCIATED(tile)) &
        & CALL finish(TRIM(routine), 'Top tile not set')

      DO WHILE (ASSOCIATED(tile))

        IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
          tile => model%Goto_next_tile()
          CYCLE
        END IF

        CALL tile%Print()

        tile => model%Goto_next_tile()

      END DO
    END IF

  END SUBROUTINE init_model_common

#endif
END MODULE mo_jsb_model_init
