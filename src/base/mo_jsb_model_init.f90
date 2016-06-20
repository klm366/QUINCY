!>
!! @brief Contains methods for JSBACH initializations.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_model_init
#ifndef __NO_JSBACH__

  USE mo_jsb_base,         ONLY: jsbach
  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: message, finish

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

  SUBROUTINE init_model_standalone(model_id, model_name)

    USE mo_jsb_grid_class,  ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,        ONLY: Register_grid
    USE mo_jsb_model_class, ONLY: print_model

    INTEGER,                    INTENT(in) :: model_id
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: model_name
!!$#ifdef HAVE_F2003
!!$    TYPE(t_jsb_model), POINTER, INTENT(inout) :: self
!!$#else
!!$    TYPE(t_jsb_model), POINTER                :: self
!!$#endif

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_standalone'

    model => jsbach%models(model_id)%m
    IF (PRESENT(model_name)) THEN
      IF (ADJUSTL(TRIM(model_name)) /= ADJUSTL(TRIM(model%name))) &
        CALL finish(TRIM(routine), 'Model ID and name are inconsistent with namelist config.')
    END IF

    ! TBD: loop through already existing grids and use it if already exists with same name
    grid => new_grid(model%config%grid_filename, model%config%grid_name)
    CALL register_grid(grid)
    model%grid_id = grid%id

    CALL print_model(model)

    NULLIFY(model)

  END SUBROUTINE init_model_standalone

!!$  SUBROUTINE init_model_patch( model_id, patch, &
!!$                             & nsoil_w, dzsoil_w, nsoil_e, dzsoil_e)
  SUBROUTINE init_model_patch( model_id, patch)

    USE mo_jsb_domain_iface, ONLY: t_patch
    USE mo_jsb_grid_class,   ONLY: t_jsb_grid, new_grid
    USE mo_jsb_grid,         ONLY: Register_grid

    INTEGER,                    INTENT(in) :: model_id
    TYPE(t_patch),              INTENT(in) :: patch

!!$    ! The soil layer structure has be sent back to ICON because this information is needed to setup the CDI output
!!$    INTEGER,           OPTIONAL, INTENT(out) :: &
!!$      & nsoil_w,                                & !< Number of soil layers (water)
!!$      & nsoil_e                                   !< Number of soil layers (energy)
!!$    REAL(wp), POINTER, OPTIONAL              :: &
!!$      & dzsoil_w(:),                            & !< Soil layer depths (water)
!!$      & dzsoil_e(:)                               !< Soil layer depths (energy)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_patch'

    IF (patch%id < 0) &
      CALL finish(TRIM(routine), 'Improper call with invalid patch')

    model => jsbach%models(model_id)%m

    grid => new_grid(patch=patch)
    CALL Register_grid(grid)
    model%grid_id = grid%id

    CALL init_model_common(model)

!!$    IF (PRESENT(nsoil_w)) nsoil_w = model%config%soil%nsoil_water
!!$    IF (PRESENT(nsoil_e)) nsoil_e = model%config%soil%nsoil_energy
!!$    IF (PRESENT(dzsoil_w)) THEN
!!$      ALLOCATE(dzsoil_w(nsoil_w))
!!$      dzsoil_w(1:nsoil_w) = model%config%soil%dz_water(1:nsoil_w)
!!$    END IF
!!$    IF (PRESENT(dzsoil_e)) THEN
!!$      ALLOCATE(dzsoil_e(nsoil_e))
!!$      dzsoil_w(1:nsoil_e) = model%config%soil%dz_energy(1:nsoil_e)
!!$    END IF

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

    model => jsbach%models(model_id)%m

    grid => new_grid(decomp=decomposition)
    CALL Register_grid(grid)
    model%grid_id = grid%id

    CALL init_model_common(model)

    NULLIFY(model)

  END SUBROUTINE init_model_decomposition

  SUBROUTINE init_model_common(model)

    USE mo_jsb_process_class,   ONLY: SEB, TURB, SOIL, HYDRO, VEG, RAD, HD, ASSIMI
    USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, ON_SUBTREE_
    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process
    USE mo_jsb_model_class,     ONLY: print_model
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_lctlib_class,    ONLY: Read_lctlib
    USE mo_jsb_model_usecases,  ONLY: init_usecase
    USE mo_seb_init,            ONLY: seb_init
    USE mo_turb_init,           ONLY: turb_init
    USE mo_soil_init,           ONLY: soil_init
    USE mo_hydro_init,          ONLY: hydro_init
    USE mo_veg_init,            ONLY: veg_init_ic, veg_init_bc
    USE mo_rad_init,            ONLY: rad_init
    ! USE mo_assimi_init,         ONLY: assimi_init
#ifndef __NO_JSBACH_HD__
    USE mo_hd_init,             ONLY: hd_init_ic, hd_init_bc
#endif
    USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,            ONLY: Register_vgrid
    USE mo_jsb_io,              ONLY: ZAXIS_SURFACE, cdiDefMissval, missval
    USE mo_util,                ONLY: one_of

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile_abstract), POINTER :: tile
    TYPE(t_jsb_vgrid),          POINTER :: surface                      ! Vertical grid
    INTEGER                             :: iproc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_model_common'

    CALL cdiDefMissval(missval)

    ALLOCATE(model%processes(max_no_of_processes))
    DO iproc=1,max_no_of_processes
      model%processes(iproc)%p => Create_process(iproc, model%id, TRIM(model%namelist_filename))
    END DO

    CALL model%Configure_processes()

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

    ! Temporary for testing ... should be in separate subroutine
    CALL model%Queue_task(RAD,    "surface_radiation")
    !CALL model%Queue_task(HYDRO, "snow_and_skin_fraction")
    CALL model%Queue_task(SEB,    "surface_energy")
    CALL model%Queue_task(SEB,    "surface_fluxes")
    CALL model%Queue_task(RAD,    "albedo")
    CALL model%Queue_task(TURB,   "surface_roughness")
    CALL model%Queue_task(ASSIMI, "assimilation_scaling_factors") ! normally called before assimilation_waterlimited
                                                                  ! and assimilation_waterunlimited
    CALL model%Queue_task(ASSIMI, "assimilation_waterlimited")    ! first call without water limitation (waterLimitationFlag =.FALSE.)
    ! R: Das ist eine grundsätzliche Entscheidung hier:
    ! R: Ein Task könnte unterschiedliche Konfigurationen habe.
    ! R: So wie die Struktur jetzt ist, MUSS jede dieser Konfigurationen
    ! R: ein eigener Task sein, denn die Queue (model%queue%next%...%task)
    ! R: zeigt immer auf einen Task der unter model%processes(2)%p%task_list(:)
    ! R: steht. Würde man zwischen den Aufrufen hier versuchen den Task zu
    ! R: konfigurieren z.B. model%queue%next%...%task%waterLimitationFlag = .TRUE.
    ! R: so würde später wenn die Queue abgearbeitet wird immer nur die letzte
    ! R: Konfiguration existieren, da der Task nur einmal existiert und die Queue
    ! R: nur Pointer auf diesen einen Task hat.
    ! R: Die Konfiguration gehört nicht zu model%processes(2)%p%config da sie
    ! R: sich nicht auf einen gesamten Prozess sondern nur auf einen Task darin
    ! R: bezieht. Außerdem wäre auch jede andere Speicherstelle die sich nicht
    ! R: bei/während der Abarbeitung der Queue ändert unnütz. Die Speicherstelle
    ! R: müsste unter model%queue% z.B. als task_config liegen.
    ! R: So wie es jetzt ist (dass die Tasks je nach Konfiguration eigenständige
    ! R: Tasks sein müssen) finde ich gar nicht so schlecht. Denn dadurch gibt
    ! R: es keine "Abhängigkeiten" des Tasks in der Queue sondern nur eindeutig
    ! R: definierte Tasks. Mit welchen IFs man die Queue auch erstellt hat, am
    ! R: Ende ist sie eindeutig und man kann sich die Tasks eindeutig ausgegeben
    ! R: lassen. Außerdem: die Konfiguration beinhaltet nicht die
    ! R: Fallunterscheidungen für jedes Tile und LCT-Typ usw....

    ! include here the tasks soil_hydrology and vegetation_conductance
    CALL model%Queue_task(ASSIMI, "assimilation_waterunlimited")  ! second call with water limitation (waterLimitationFlag =.TRUE.)
                                                   ! should be called after the soil hydrology is calculated (former update_soil_1)
                                                   !              and after the vegetation conductance is calculated
                                                         !              (former update_vegetation_conductanc)
    CALL model%Queue_task(TURB,   "humidity_scaling")
    CALL model%Queue_task(SEB,    "asselin_filter")      ! Should always be the last task

    CALL print_model(model)
    CALL model%Print_actions()

    ! At this point, the model configuration, tile and memory structure is set up.

    ! We continue with initializing the different model processes.

    ! Iterate over tiles and call initialization if necessary
    CALL model%Reset_tiles()
    tile => model%Get_top_tile()
    IF (.NOT. ASSOCIATED(tile)) &
      & CALL finish(TRIM(routine), 'Top tile not set')

    DO WHILE (ASSOCIATED(tile))

      IF (tile%visited) THEN  ! Have been here before
        tile => model%Goto_next_tile()
        CYCLE
      END IF

      ! Initialize tile memory and read tile fractions
!!$      IF (ASSOCIATED(tile%parent_tile)) THEN
        CALL tile%Init(TRIM(model%shortname), '', TRIM(tile%name))
!!$      ELSE
!!$        CALL tile%Init(TRIM(model%shortname), '', '')
!!$      END IF

      ! TBD: Loop over processes instead of explicit
      IF (one_of(tile%process_action(RAD), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        CALL rad_init(tile)
      END IF

      IF (one_of(tile%process_action(SEB), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        CALL seb_init(tile)
      END IF

      IF (one_of(tile%process_action(TURB), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        CALL turb_init(tile)
      END IF

      IF (one_of(tile%process_action(SOIL), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        CALL soil_init(tile)
      END IF

      IF (one_of(tile%process_action(HYDRO), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        CALL hydro_init(tile)
      END IF

      ! Vegetation must be initialized after soil!
      IF (one_of(tile%process_action(VEG), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        CALL veg_init_bc(tile)
        CALL veg_init_ic(tile)
      END IF

      ! IF (one_of(tile%process_action(ASSIMI), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
      !   CALL assimi_init_bc(tile)
      !   CALL assimi_init_ic(tile)
      ! END IF

#ifndef __NO_JSBACH_HD__
      IF (one_of(tile%process_action(HD), (/ON_LEAFS_,ON_TILE_,ON_SUBTREE_/)) > 0) THEN
        CALL hd_init_bc(tile)
        CALL hd_init_ic(tile)
      END IF
#endif

      tile => model%Goto_next_tile()

    END DO

    ! Initialize tile fractions for all tiles
    ! This is done after the initialization of the initial and boundary conditions of the model components
    ! because e.g. the fraction of the vegetated tile depends on LAI.
    ! TBD: later read veg_ratio_max into 'veg' tile, and set 'bare' tile to 1-veg_ratio_max

  END SUBROUTINE init_model_common

#endif
END MODULE mo_jsb_model_init
