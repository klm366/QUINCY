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
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_model_init
#ifndef __NO_JSBACH__

  dsl4jsb_Use_processes SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_, CARBON_, DISTURB_, FUEL_
  dsl4jsb_Use_processes ALCC_, A2L_, VEG_, DYNVEG_
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_config(DYNVEG_)

  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_control,      ONLY: get_no_of_models, debug_on
  USE mo_jsb_model_class,  ONLY: t_jsb_model

  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: finish
  USE mo_jsb_parallel,     ONLY: Get_omp_thread, Is_omp_inside_serial
  USE mo_jsb_time,         ONLY: is_time_restart, get_time_start

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_setup_grid, jsbach_init, jsbach_init_after_restart

  INTERFACE jsbach_setup_grid
    MODULE PROCEDURE setup_grid_standalone
    MODULE PROCEDURE setup_grid_patch
  END INTERFACE jsbach_setup_grid

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_init'

CONTAINS

  SUBROUTINE setup_grid_standalone(nproca, nprocb, npedim, model_id, model_name)

#ifndef __ICON__
    !>>> TODO: echam
    USE mo_gaussgrid,     ONLY : inigau
    USE mo_echam_convect_tables, ONLY : init_convect_tables
    USE mo_netCDF,      ONLY: add_dim
    !<<<
    USE mo_jsb_grid_iface,   ONLY: get_lons_and_lats_from_file, get_nlon_nlat_from_file
    USE mo_jsb_domain_iface, ONLY: t_patch, jsbalone_init_decomp
    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,         ONLY: Register_grid
#endif

    INTEGER,                    INTENT(in) :: model_id, nproca, nprocb, npedim
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: model_name

#ifndef __ICON__
    TYPE(t_jsb_model),  POINTER  :: model
    TYPE(t_jsb_grid),   POINTER  :: grid
    TYPE(t_patch)                :: decomposition

    INTEGER :: nlon, nlat
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':setup_grid_standalone'

#ifdef __ICON__
    ! avoid compiler warnings about dummy arguments not being used
    IF (.FALSE.) WRITE(*,*) nproca, nprocb, npedim, model_id, model_name

    CALL finish(TRIM(routine),'This subroutine should not be called with ICON')
#else
    model => Get_model(model_id)

    IF (PRESENT(model_name)) THEN
      IF (ADJUSTL(TRIM(model_name)) /= ADJUSTL(TRIM(model%name))) &
        CALL finish(TRIM(routine), 'Model ID and name are inconsistent with namelist config.')
    END IF

    CALL get_nlon_nlat_from_file(model%config%grid_filename, nlon, nlat)
    CALL jsbalone_init_decomp(nlon, nlat, nproca, nprocb, npedim, decomposition)

    CALL inigau               ! beside others required for get_area, i.e. required before new_grid
    CALL init_convect_tables

    CALL add_dim("lon", nlon, "longitude", "degrees_east" ) !required to write restart files
    CALL add_dim("lat", nlat, "latitude", "degrees_north")

    !needs to be done after the decomposition is initialised; required for get_lon and get_lat, especially for site-level sims
    CALL get_lons_and_lats_from_file(model%config%grid_filename, decomposition%nproma, decomposition%ngpblks)
    grid => new_grid(decomposition, type='echam')

    CALL register_grid(grid)
    model%grid_id = grid%id

    ! CALL init_model_common(model)

    NULLIFY(model)
#endif

  END SUBROUTINE setup_grid_standalone

  SUBROUTINE setup_grid_patch( model_id, patch, type)

    USE mo_jsb_domain_iface, ONLY: t_patch
    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,         ONLY: Register_grid

    INTEGER,               INTENT(in) :: model_id
    TYPE(t_patch), TARGET, INTENT(in) :: patch
    CHARACTER(LEN=*),      INTENT(in) :: type

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':setup_grid_patch'

    IF (model_id > get_no_of_models()) &
      CALL finish(routine, 'Host request to more jsbach models than are defined')

    model => Get_model(model_id)

    grid => new_grid(patch, type)
    CALL Register_grid(grid)
    model%grid_id = grid%id

    ! CALL init_model_common(model)

    NULLIFY(model)

  END SUBROUTINE setup_grid_patch

  SUBROUTINE jsbach_init(model_id)

    USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_, AGGREGATE_
    USE mo_jsb_model_class,     ONLY: print_model
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_lct_class,       ONLY: Contains_lct
    USE mo_seb_init,            ONLY: seb_init
    USE mo_turb_init,           ONLY: turb_init
    USE mo_sse_init,            ONLY: sse_init
    USE mo_hydro_init,          ONLY: hydro_init
    USE mo_rad_init,            ONLY: rad_init
    USE mo_atm2land_init,       ONLY: atm2land_init
    USE mo_assimi_init,         ONLY: assimi_init
    USE mo_pheno_init,          ONLY: pheno_init
    USE mo_veg_init,            ONLY: veg_init
    USE mo_carbon_init,         ONLY: carbon_init
    USE mo_fuel_init,           ONLY: fuel_init
    !USE mo_disturb_init,        ONLY: disturb_init
#ifndef __NO_JSBACH_HD__
    USE mo_hd_init,             ONLY: hd_init_ic, hd_init_bc
#endif
    USE mo_jsb_io,              ONLY: cdiDefMissval, missval
    USE mo_util,                ONLY: one_of

    INTEGER, INTENT(in) :: model_id

    TYPE(t_jsb_model), POINTER  :: model
    CLASS(t_jsb_tile_abstract), POINTER :: tile, child
    INTEGER                             :: ilct, jlct, i, no_children
    REAL(wp), ALLOCATABLE               :: fract_tmp(:,:)
    INTEGER                             :: no_omp_thread

    dsl4jsb_Def_config(PHENO_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_init'

    no_omp_thread = Get_omp_thread()

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    CALL cdiDefMissval(missval)

    model => Get_model(model_id)

    dsl4jsb_Get_config(PHENO_)

    !!---------------------------------------------------------------------------------------------
    !! Define Task Queue 
    !!   (Temporary for testing ... should be in separate subroutine !!)
    !! 
    !! Two Queues: default jsb4 & quincy
    IF (.NOT. model%config%use_quincy) THEN
      !NOTE: JN changed the order of the processes
      !- processes potentially changing the area of a tile are conducted first in the timestep,
      !- but the task npp_buffer is performed even before to provide already updated multi annual average NPP to DYNVEG at the beginning of a year
      !- PHENO processes were moved up: need to happen before C_NPP_pot_allocation for allometric considerations
      CALL model%Queue_task(ASSIMI_, "npp_buffer")
      !!------------------------ tasks that can change area
      CALL model%Queue_task(DYNVEG_,"dynveg")
      !! 
      !Disturbances happen on a daily basis and can lead to area changes when running with cohorts
      IF (model%processes(FUEL_)%p%config%active) THEN
        CALL model%Queue_task(FUEL_, "fuel")
      END IF
      IF (model%processes(DISTURB_)%p%config%active) THEN
         CALL model%Queue_task(DISTURB_, "natural_disturbances")
      END IF
  
      !!------------------------ tasks that do not change area
      !!
      CALL model%Queue_task(RAD_,    "surface_radiation")
      CALL model%Queue_task(PHENO_,  "phenology")           ! JN moved phenology up - for allom it needs to happen
                                                            ! before C_NPP_pot_allocation
                                                            ! and with area changes it should happen before assimi
                        ! --- before it was between humidity_scaling and foliage_projected_cover
                        ! + prev comment: can be placed between the task "surface_radiation" and "asselin_filter"
      CALL model%Queue_task(PHENO_,  "foliage_projected_cover")
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
      CALL model%Queue_task(SSE_,    "soil_and_snow_properties")
      CALL model%Queue_task(SSE_,    "soil_and_snow_temperature")
      CALL model%Queue_task(HYDRO_,  "soil_hydrology")
#ifndef __NO_JSBACH_HD__
      CALL model%Queue_task(HD_,     "hd_local_outflow")
#endif
      CALL model%Queue_task(SEB_,    "snowmelt_correction")
      IF (model%processes(ASSIMI_)%p%config%active) THEN
        ! Normally called before update_canopy_cond_unstressed and update_assimilation
        CALL model%Queue_task(ASSIMI_, "assimilation_scaling_factors")
      END IF
      IF (model%processes(ASSIMI_)%p%config%active) THEN
        CALL model%Queue_task(ASSIMI_, "canopy_cond_unstressed")
      ELSE
        CALL model%Queue_task(HYDRO_,  "canopy_cond_unstressed")
      END IF
      CALL model%Queue_task(HYDRO_, "water_stress")
      IF (model%processes(ASSIMI_)%p%config%active) THEN
        CALL model%Queue_task(ASSIMI_, "canopy_cond_stressed")
        CALL model%Queue_task(ASSIMI_, "assimilation")
      ELSE
        CALL model%Queue_task(HYDRO_,  "canopy_cond_stressed")
      END IF
      ! CALL model%Queue_task(HYDRO_,  "snow_and_skin_fraction")
      IF (model%processes(CARBON_)%p%config%active) THEN
         CALL model%Queue_task(CARBON_, "C_NPP_pot_allocation")
      END IF
      CALL model%Queue_task(TURB_,   "humidity_scaling")
      CALL model%Queue_task(TURB_,   "surface_roughness")
      ! Note: the albedo is updated at the end of the time step, one time step before the radiation is updated in the
      !       atmosphere. Otherwise, in between, the albedo changes would be consindered in the jsbach radiation but where
      !       not considered in the atmospheric radiation calculations. This would be inconsistent. Furthermore, atmospheric
      !       cloud changes affect radiation also only at radiation time step and jsbach should be also consistent with this.
      !       In the albedo update time step the albedo should be updated after the SEB process, sothat also the SEB remains
      !       consistent within this time step.
      CALL model%Queue_task(RAD_,    "albedo")
      ! Note: The explanation above for the aldedo accounts also for alcc as it may also changes surface properties
      !       affecting atmospheric radiation. However, in alcc we do not redistribute the albedos of the pfts because
      !       they shall not keep their quantity for the gridbox but change with alcc. Therefore, the sequence between the alcc and
      !       albedo process does not matter.
      IF (model%processes(ALCC_)%p%config%active) THEN
        CALL model%Queue_task(ALCC_, "alcc")
      END IF
      CALL model%Queue_task(HYDRO_,  "water_balance")
      CALL model%Queue_task(SEB_,    "asselin_filter")      ! Should always be the last task
    
    !!---------------------------------------------------------------------------------------------
    !! quincy
    !!  usecase must be set to: "jsbach_pfts" (which assumes ASSIMI_ to be active)      (this is not tested)
    !!  CARBON_ may be active                                                           (this is not tested)
    !!  l_forestRegrowth is assumed to be .FALSE.                                       (this is not tested)
    ELSE IF (model%config%use_quincy) THEN
      !!------------------------ tasks that simply set fluxes back to zero at the beginning of each timestep
      !!
      CALL model%Queue_task(ASSIMI_, "reset_assimi_fluxes")         !! quincy
      CALL model%Queue_task(VEG_,    "reset_veg_fluxes")            !! quincy
      
      !!------------------------ tasks that can change area (e.g. disturbance & fuel)
      !! 
      ! !Disturbances happen on a daily basis and can lead to area changes when running with cohorts
      ! IF (model%processes(FUEL_)%p%config%active) THEN
      !   CALL model%Queue_task(FUEL_, "fuel")
      ! END IF
      ! IF (model%processes(DISTURB_)%p%config%active) THEN
      !    CALL model%Queue_task(DISTURB_, "natural_disturbances")
      ! END IF
  
      !!------------------------ tasks that do not change area
      !!
      CALL model%Queue_task(RAD_,    "sw_radiation_q_")             !! quincy (sw radiation and albedo)
      CALL model%Queue_task(RAD_,    "surface_radiation")
      CALL model%Queue_task(PHENO_,  "phenology_q_")                !! quincy
      CALL model%Queue_task(PHENO_,  "foliage_projected_cover")
      CALL model%Queue_task(HYDRO_,  "snow_and_skin_fraction")      ! may be called two times
      CALL model%Queue_task(HYDRO_,  "soil_properties")
      CALL model%Queue_task(SEB_,    "surface_energy")
      CALL model%Queue_task(HYDRO_,  "evaporation")
      CALL model%Queue_task(SEB_,    "surface_fluxes")
      CALL model%Queue_task(HYDRO_,  "surface_hydrology")
      CALL model%Queue_task(HYDRO_,  "snow_and_skin_fraction")      ! may be called two times
      CALL model%Queue_task(SSE_,    "soil_and_snow_properties")
      CALL model%Queue_task(SSE_,    "soil_and_snow_temperature")
      CALL model%Queue_task(HYDRO_,  "soil_hydrology")
      CALL model%Queue_task(SEB_,    "snowmelt_correction")
      CALL model%Queue_task(HYDRO_,  "water_stress")
      CALL model%Queue_task(ASSIMI_, "update_canopy_fluxes")        !! quincy
      CALL model%Queue_task(HYDRO_,  "snow_and_skin_fraction")
      CALL model%Queue_task(TURB_,   "humidity_scaling")
      CALL model%Queue_task(TURB_,   "surface_roughness")
      ! Note: albedo is updated at the end of the time step, one time step before the radiation is updated in the atmosphere
      CALL model%Queue_task(RAD_,    "albedo")
      CALL model%Queue_task(HYDRO_,  "water_balance")
      CALL model%Queue_task(SEB_,    "asselin_filter")              ! Should always be the last task
      
      ! averaging variables in time at the end of the timestep
      CALL model%Queue_task(PHENO_,  "time_average_phenology")      !! quincy
    ENDIF
    !! Task Queue created
    !!---------------------------------------------------------------------------------------------    
    
    ! debug message after creating task queue
    IF (debug_on()) THEN
      CALL print_model(model)
      CALL model%Print_actions()
    END IF

    ! At this point, the model configuration and tile structure is set up.

    ! We continue with initializing the memory structure and different model processes.

    ! Iterate over tiles and call initialization if necessary
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)

    IF (.NOT. ASSOCIATED(tile)) CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      ! Initialize tile process memories and read tile fractions
      ! @todo: once we have land use change or dynveg, l_fixed_fractions needs to be adapted
      CALL tile%Init(TRIM(model%shortname), '', TRIM(tile%name),                &
        & grid_id = model%grid_id,                                              &
        & in_var_groups=one_of(TRIM(tile%name), model%config%output_tiles) > 0)

      CALL tile%Init_fractions(                                           &
        & TRIM(model%shortname), '', TRIM(tile%name),                     &
        & l_fixed_fractions=.NOT. model%processes(ALCC_)%p%config%active, &
        & l_rel_fractions=model%config%relative_fractions_in_file)

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

      ! Initialize PHENO_ before HYDRO_ and VEG_
      IF (one_of(tile%process_action(PHENO_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL pheno_init(tile)
      END IF

      IF (one_of(tile%process_action(VEG_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL veg_init(tile)
      END IF

      ! HYDRO_ needs to be initialized before SSE_
      IF (one_of(tile%process_action(HYDRO_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL hydro_init(tile)
      END IF

      IF (one_of(tile%process_action(SSE_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL sse_init(tile)
      END IF

      IF (one_of(tile%process_action(A2L_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL atm2land_init(tile)
      END IF
      
      IF (one_of(tile%process_action(ASSIMI_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL assimi_init(tile)
        ! CALL assimi_init_bc(tile)
        ! CALL assimi_init_ic(tile)
      END IF

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

#ifndef __NO_JSBACH_HD__
      IF (one_of(tile%process_action(HD_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        CALL hd_init_bc(tile)
        CALL hd_init_ic(tile)
      END IF
#endif

      CALL model%Goto_next_tile(tile)
    END DO

    ! Iterate over tiles again and set lct fractions (after ALL tile fractions have been set above)
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)

    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    ALLOCATE(fract_tmp(SIZE(tile%lcts(1)%fract,1), SIZE(tile%lcts(1)%fract,2)))

    DO WHILE (ASSOCIATED(tile))

      IF (.NOT. tile%visited(no_omp_thread) .AND. tile%Has_children()) THEN
        CALL model%Goto_next_tile(tile)
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
              CALL tile%Get_fraction(fract=fract_tmp(:,:))
              tile%lcts(ilct)%fract(:,:) = MERGE(tile%lcts(ilct)%fract(:,:) / fract_tmp(:,:), 0._wp, fract_tmp(:,:) > 0._wp)
            END IF
            child => child%Get_next_sibling_tile()
          END DO
        ELSE
          IF (ilct > 1) CALL finish(TRIM(routine), 'Leave tiles should only have one LCT')
          tile%lcts(1)%fract(:,:) = 1._wp        ! Leaves have only one LCT
        END IF
        !$acc update device(tile%lcts(ilct)%fract)
      END DO

      CALL model%Goto_next_tile(tile)

    END DO

    DEALLOCATE(fract_tmp)

    ! Iterate over tiles again and initialize hierarchical structure for memories and variables
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (.NOT. tile%visited(no_omp_thread) .AND. tile%Has_children()) THEN
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      CALL tile%Init_vars()
      CALL tile%Cache_GPU_pointers()

      CALL model%Goto_next_tile(tile)

    END DO

    ! Iterate over tiles and crossckeck tile structure
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      IF (one_of(tile%process_action(PHENO_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        dsl4jsb_Get_config(PHENO_)
        IF (      tile%is_vegetation &
          & .AND. TRIM(dsl4jsb_Config(PHENO_)%scheme) /= 'climatology' &
          & .AND. tile%lcts(1)%lib_id == 0) &
            & CALL finish(TRIM(routine), 'Must use phenology climatology on general, non-PFT VEG tile '//TRIM(tile%name))
      END IF

      CALL model%Goto_next_tile(tile)

    END DO
    !
    !
    IF (debug_on()) THEN

      ! Iterate over tiles again and print tile information (after ALL tiles have been initialized above)
      CALL model%Reset_tiles()

      CALL model%Get_top_tile(tile)

      IF (.NOT. ASSOCIATED(tile)) &
        & CALL finish(TRIM(routine), 'Top tile not set')

      DO WHILE (ASSOCIATED(tile))

        IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
          CALL model%Goto_next_tile(tile)
          CYCLE
        END IF

        CALL tile%Print()

        CALL model%Goto_next_tile(tile)

      END DO
    END IF

!#ifdef _OPENACC
!    IF (.NOT. is_time_restart(get_time_start())) CALL model%cpu_to_gpu()
!#endif

  END SUBROUTINE jsbach_init

  SUBROUTINE jsbach_init_after_restart(model_id)

    USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_, AGGREGATE_
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_carbon_init,         ONLY: carbon_read_cpools
    USE mo_util,                ONLY: one_of

    INTEGER, INTENT(in) :: model_id

    TYPE(t_jsb_model), POINTER  :: model
    CLASS(t_jsb_tile_abstract), POINTER :: tile
    INTEGER                             :: no_omp_thread

    dsl4jsb_Def_config(CARBON_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_init_after_restart'

    no_omp_thread = Get_omp_thread()

    IF (.NOT. Is_omp_inside_serial()) THEN
      CALL finish(TRIM(routine), 'Should not be called within parallel OMP region')
    END IF

    model => Get_model(model_id)

    dsl4jsb_Get_config(CARBON_)

    ! Iterate over tiles and call initialization if necessary
    CALL model%Reset_tiles()
    CALL model%Get_top_tile(tile)
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        call model%Goto_next_tile(tile)
        CYCLE
      END IF

      ! TBD: Loop over processes instead of explicit
      IF (one_of(tile%process_action(CARBON_), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_,AGGREGATE_/)) > 0) THEN
        IF (dsl4jsb_Config(CARBON_)%read_cpools) CALL carbon_read_cpools(tile)
      END IF

      CALL model%Goto_next_tile(tile)

    END DO

#ifdef _OPENACC
    CALL model%cpu_to_gpu()
#endif

  END SUBROUTINE jsbach_init_after_restart

#endif
END MODULE mo_jsb_model_init
