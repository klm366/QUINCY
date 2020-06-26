!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains the memory class for the soil_biogeochemistry processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20170103)
!!
MODULE mo_sb_memory_class

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable
  USE mo_jsb_impl_constants,     ONLY: false
  
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE
  

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_memory, max_no_of_vars
  PUBLIC :: t_sb_main_pool, t_sb_tissues_pool, t_sb_elements_pool

  INTEGER, PARAMETER :: max_no_of_vars      = 480  !> TYPE t_jsb_var | 'normal' var: ~ 200 | 4 pools: 240 var | 2 tissues: 12

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb_elements_pool
  !! pool with only variables
  !! last pool in line
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_pool) :: t_sb_elements_pool
    TYPE(t_jsb_var_real3d), POINTER :: carbon
    TYPE(t_jsb_var_real3d), POINTER :: nitrogen
    TYPE(t_jsb_var_real3d), POINTER :: phosphorus
    TYPE(t_jsb_var_real3d), POINTER :: carbon13
    TYPE(t_jsb_var_real3d), POINTER :: carbon14
    TYPE(t_jsb_var_real3d), POINTER :: nitrogen15
  CONTAINS
    PROCEDURE :: Init => t_sb_elements_pool_init
    PROCEDURE :: Get_element_name_by_id => t_sb_elements_pool_get_element_name_by_id
  END TYPE t_sb_elements_pool

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for t_sb_tissues_pool
  !! vegetation pool of tissues
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_sb_elements_pool) :: t_sb_tissues_pool
    TYPE(t_sb_elements_pool), POINTER :: dom                 !< elements of dissolved organic matter
    TYPE(t_sb_elements_pool), POINTER :: dom_assoc           !< elements of minerally associated disolved organic matter
    TYPE(t_sb_elements_pool), POINTER :: soluable_litter     !< elements of metabolic litter (sugars, cellulose etc.)
    TYPE(t_sb_elements_pool), POINTER :: polymeric_litter    !< elements of structural litter (lignified litter)
    TYPE(t_sb_elements_pool), POINTER :: woody_litter        !< elements of physically protected wood litter
    TYPE(t_sb_elements_pool), POINTER :: fungi               !< elements of saprophytic fungal biomass
    TYPE(t_sb_elements_pool), POINTER :: mycorrhiza          !< elements of mycorrhizal fungal biomass 
    TYPE(t_sb_elements_pool), POINTER :: microbial           !< elements of microbial biomass (bacteria!)
    TYPE(t_sb_elements_pool), POINTER :: residue             !< elements of necromass of fungi, microbial and mycorrhiza
    TYPE(t_sb_elements_pool), POINTER :: residue_assoc       !< elements of minerally associated microbial necromass 
  END TYPE t_sb_tissues_pool

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for t_sb_main_pool
  !! pool-structure containing pools and fluxes
  !! 1st pool in line
  !! added directly into the t_sb_memory as "CLASS(t_jsb_pool), POINTER :: pools"
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_sb_elements_pool) :: t_sb_main_pool
    ! pools
    TYPE(t_sb_tissues_pool),   POINTER :: sb_pool               !< 
    TYPE(t_sb_elements_pool),  POINTER :: sb_delta_pools        !< 
    ! fluxes
    TYPE(t_sb_tissues_pool),   POINTER :: sb_formation          !< 
    TYPE(t_sb_tissues_pool),   POINTER :: sb_loss               !< 
    TYPE(t_sb_tissues_pool),   POINTER :: sb_transport          !< 
    TYPE(t_sb_elements_pool),  POINTER :: sb_mycorrhiza_export  !< 
  END TYPE t_sb_main_pool

  !-----------------------------------------------------------------------------------------------------
  !> list with the element_id of all elements as defined above
  !! is allocated in Init_PROC_memory
  !-----------------------------------------------------------------------------------------------------
  INTEGER, ALLOCATABLE, SAVE :: all_elements(:)

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb memory
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_sb_memory

    !> pools & fluxes
    !! 
    ! the sb_main_pool does contain all other pools & fluxes 
    TYPE(t_sb_main_pool), POINTER :: sb_main_pool

    ! 2D total fluxes for mass balance checks in soil_biogeochemistry model
    TYPE(t_jsb_var_real2d)          :: & 
                                    total_flux_carbon, &          !< carbon mass balance checks in soil_biogeochemistry model [mol C m-2 s-1]
                                    total_flux_nitrogen, &        !< nitrogen mass balance checks in soil_biogeochemistry model [mol N m-2 s-1]
                                    total_flux_phosphorus, &      !< phosphorus mass balance checks in soil_biogeochemistry model [mol P m-2 s-1]
                                    total_flux_carbon13, &        !< carbon13 mass balance checks in soil_biogeochemistry model [mol C13 m-2 s-1]
                                    total_flux_carbon14, &        !< carbon14 mass balance checks in soil_biogeochemistry model [mol C14 m-2 s-1]
                                    total_flux_nitrogen15         !< nitrogen15 mass balance checks in soil_biogeochemistry model [mol N15 m-2 s-1]

    ! 3D soil properties
    TYPE(t_jsb_var_real3d)          :: &
                                    ph_sl                         !< ph per soil layer

  CONTAINS
    PROCEDURE :: Init => Init_sb_memory
  END TYPE t_sb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> helper function for t_sb_elements_pool
  !! PROCEDURE :: Get_element_name_by_id => t_sb_elements_pool_get_element_name_by_id
  !-----------------------------------------------------------------------------------------------------
  FUNCTION t_sb_elements_pool_get_element_name_by_id(this, id) RESULT(name)

    CLASS(t_sb_elements_pool),  INTENT(in)  :: this
    INTEGER,                    INTENT(in)  :: id
    CHARACTER(LEN=:),           ALLOCATABLE :: name

    ! name may be consistent with the "name" provided to "CALL this%Add_element()"
    SELECT CASE (id)
    CASE (ELEM_C)
      name='carbon'
    CASE (ELEM_N)
      name='nitrogen'
    CASE (ELEM_P)
      name='phosphorus'
    CASE (ELEM_C13)
      name='carbon13'
    CASE (ELEM_C14)
      name='carbon14'
    CASE (ELEM_N15)
      name='nitrogen15'
    END SELECT

  END FUNCTION t_sb_elements_pool_get_element_name_by_id

  !-----------------------------------------------------------------------------------------------------
  !> initializing an object of the SB_ elements pool "t_sb_elements_pool"
  !! "pool_object%Init()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_elements_pool_init(this, name, shortname, l_elements, element_list, element_unit)

    CLASS(t_sb_elements_pool),  INTENT(inout) :: this               ! object of TYPE t_sb_elements_pool
    CHARACTER(LEN=*),           INTENT(in)    :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),           INTENT(in)    :: shortname          ! pool short name (No spaces !)
    LOGICAL,          OPTIONAL, INTENT(in)    :: l_elements         ! IF pool has elements 
    INTEGER,          OPTIONAL, INTENT(in)    :: element_list(:)    ! list of elements; use if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: element_unit       ! unit of elements

    LOGICAL                     :: l_add_elements
    INTEGER, ALLOCATABLE        :: element_list_loc(:)
    INTEGER                     :: i
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_elements_pool_init'

    !--------------------------------------------
    this%name = name
    this%shortname = shortname

    l_add_elements = .TRUE.
    IF (PRESENT(l_elements)) l_add_elements = l_elements
    IF (.NOT. l_add_elements) RETURN

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')

    ALLOCATE(element_list_loc(SIZE(all_elements)))
    element_list_loc(:) = all_elements(:)
    CALL assign_if_present_allocatable(element_list_loc, element_list)  ! selects element_list if available
  
    this%element_unit = element_unit
  
    !> carbon
    IF (One_of(ELEM_C, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=3)
      i = this%Find_element_position(ELEM_C)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon => selector
        END SELECT
      END IF
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=3)
      i = this%Find_element_position(ELEM_N)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%nitrogen => selector
        END SELECT
      END IF
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=3)
      i = this%Find_element_position(ELEM_P)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%phosphorus => selector
        END SELECT
      END IF
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=3)
      i = this%Find_element_position(ELEM_C13)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon13 => selector
        END SELECT
      END IF
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=3)
      i = this%Find_element_position(ELEM_C14)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon14 => selector
        END SELECT
      END IF
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=3)
      i = this%Find_element_position(ELEM_N15)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%nitrogen15 => selector
        END SELECT
      END IF
    END IF

  END SUBROUTINE t_sb_elements_pool_init

  !-----------------------------------------------------------------------------------------------------
  !> create and init an object of the SB_ tissues pool "t_sb_tissues_pool"
  !! init and add all its objects of the "t_sb_elements_pool" with "pool_object%Init()"
  !! "object => Create_sb_tissues_pool()"
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_sb_tissues_pool(name, shortname, unit, lct_ids, pool_prefix, pool_suffix) RESULT(sb_tissues_pool)

    CHARACTER(LEN=*),           INTENT(in)  :: name
    CHARACTER(LEN=*),           INTENT(in)  :: shortname
    CHARACTER(LEN=*),           INTENT(in)  :: unit
    INTEGER,                    INTENT(in)  :: lct_ids(:)
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)  :: pool_prefix
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)  :: pool_suffix
    TYPE(t_sb_tissues_pool),    POINTER     :: sb_tissues_pool

    TYPE(t_sb_elements_pool),   POINTER     :: sb_elements_pool
    CHARACTER(LEN=:),           ALLOCATABLE :: prefix_loc, suffix_loc

    prefix_loc = ''
    suffix_loc = ''
    CALL assign_if_present_allocatable(prefix_loc, pool_prefix)
    CALL assign_if_present_allocatable(suffix_loc, pool_suffix)
    IF (prefix_loc /= '') prefix_loc = prefix_loc // '_'
    IF (suffix_loc /= '') suffix_loc = '_' // suffix_loc

    ! tissue pool init
    ALLOCATE(sb_tissues_pool)
    CALL sb_tissues_pool%Init(name, shortname, l_elements=.TRUE., element_unit=unit)
    
    ! element pools allocate & init & add to tissue pool
    ! "shortname" is used for model output
    ! "element_list" is optional and can be used to select specific elements (not needed when all elements may be added)
    ! dom
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='dom_pool'//suffix_loc, shortname='dom', l_elements=.TRUE., element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%dom => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! dom_assoc
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='dom_assoc_pool'//suffix_loc, shortname='dom_assoc', l_elements=.TRUE., element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%dom_assoc => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! soluable_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='soluable_litter_pool'//suffix_loc, shortname='sol_litter', l_elements=.TRUE., &
                               element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%soluable_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! polymeric_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='polymeric_litter_pool'//suffix_loc, shortname='pol_litter', l_elements=.TRUE., &
                               element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%polymeric_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! woody_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='woody_litter_pool'//suffix_loc, shortname='woo_litter', l_elements=.TRUE., element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%woody_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! fungi
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='fungi_pool'//suffix_loc, shortname='fungi', l_elements=.TRUE., element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%fungi => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! mycorrhiza
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='mycorrhiza_pool'//suffix_loc, shortname='mycorrhiza', l_elements=.TRUE., element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%mycorrhiza => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! microbial
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='microbial_pool'//suffix_loc, shortname='microbial', l_elements=.TRUE., element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%microbial => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! residue
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='residue_pool'//suffix_loc, shortname='residue', l_elements=.TRUE., element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%residue => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! residue_assoc
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='residue_assoc_pool'//suffix_loc, shortname='residue_assoc', l_elements=.TRUE., &
                               element_unit=unit)
    CALL sb_tissues_pool%Add_pool(sb_elements_pool)
    sb_tissues_pool%residue_assoc => sb_elements_pool
    NULLIFY(sb_elements_pool)

  END FUNCTION Create_sb_tissues_pool

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SB_ process
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_sb_memory),   INTENT(inout), TARGET  :: mem             !> SB_ memory
    CHARACTER(len=*),     INTENT(in)             :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)             :: suffix          !> tile name
    INTEGER,              INTENT(in)             :: lct_ids(:)      !> vector of land cover type IDs
    INTEGER,              INTENT(in)             :: model_id        !> model ID model\%id
    ! ---------------------------
    TYPE(t_jsb_model),          POINTER :: model
    TYPE(t_jsb_grid),           POINTER :: hgrid                    ! Horizontal grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_surface            ! Vertical grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_soil_sb            ! Vertical grid
    TYPE(t_sb_elements_pool),   POINTER :: sb_elements_pool         ! pools/fluxes like sb_delta_pools & sb_mycorrhiza_export
    TYPE(t_sb_tissues_pool),    POINTER :: sb_tissues_pool          ! pools/fluxes like sb_pool & sb_formation
    INTEGER                     :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_memory'
    ! ---------------------------
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')


    !> add memory only for LAND & PFT tiles
    !!
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      !=============================================================================================================================
      !> pools & fluxes
      !=============================================================================================================================

      !-----------------------------------------------------------------------------------------------------------------------------
      !> populate all_elements(:) with the element_id of all elements
      !! QUINCY is using always all of the below six elements
      IF (.NOT. ALLOCATED(all_elements)) THEN
        ALLOCATE(all_elements(6))
        all_elements(:) = [ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15]
      END IF

      !-----------------------------------------------------------------------------------------------------------------------------
      !> create and populate the main pool structure mem%pools that contains all other pools
      !! note: no spaces in 'name' and/or 'short name'
      ALLOCATE(t_sb_main_pool :: mem%pools)
      CALL mem%pools%Init(name='soil_biogeochemistry_main_pool', shortname='sb_main_pool', l_elements=.FALSE.)

      SELECT TYPE (selected => mem%pools)
      CLASS IS (t_sb_main_pool)
        mem%sb_main_pool => selected
      END SELECT

      !-----------------------------------------------------------------------------------------------------------------------------
      !> create the vegetation tissue pools & fluxes (sb_tissues_pool) 
      !! sb_pool / sb_formation / ...
      ! create pools & fluxes
      ! "shortname" is used for model output names
      ! sb_pool
      sb_tissues_pool => Create_sb_tissues_pool(name='soil_biogeochemistry_pool', shortname='sb_pool', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbpool') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_tissues_pool)
      mem%sb_main_pool%sb_pool => sb_tissues_pool
      ! sb_formation
      sb_tissues_pool => Create_sb_tissues_pool(name='sb_formation_flux', shortname='sb_formation', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbform') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_tissues_pool)
      mem%sb_main_pool%sb_formation => sb_tissues_pool
      ! sb_loss
      sb_tissues_pool => Create_sb_tissues_pool(name='sb_loss_flux', shortname='sb_loss', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbloss') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_tissues_pool)
      mem%sb_main_pool%sb_loss => sb_tissues_pool
      ! sb_transport
      sb_tissues_pool => Create_sb_tissues_pool(name='sb_transport_flux', shortname='sb_transport', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbtran') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_tissues_pool)
      mem%sb_main_pool%sb_transport => sb_tissues_pool
      NULLIFY(sb_tissues_pool)

      ! init and add element pools to the veg_main_pool
      ! (i.e., add them at the same level as the above veg_tissues_pool veg_pool / veg_growth / ...)
      ! "shortname" is used for model output names
      ! sb_delta_pools
      ALLOCATE(sb_elements_pool)
      CALL sb_elements_pool%Init(name='sb_delta_pools_elem_pool', shortname='sb_delta_pools', &
        &                        l_elements=.TRUE., element_unit='mol m-3')
      CALL mem%sb_main_pool%Add_pool(sb_elements_pool)
      mem%sb_main_pool%sb_delta_pools => sb_elements_pool
      NULLIFY(sb_elements_pool)
      ! sb_mycorrhiza_export
      ALLOCATE(sb_elements_pool)
      CALL sb_elements_pool%Init(name='sb_mycorrhiza_export_elem_flux', shortname='sb_mycorrhiza_export', &
        &                        l_elements=.TRUE., element_unit='mol m-3')
      CALL mem%sb_main_pool%Add_pool(sb_elements_pool)
      mem%sb_main_pool%sb_mycorrhiza_export => sb_elements_pool
      NULLIFY(sb_elements_pool)

      ! set paths of all above initialized & added pools (type-bound procedure of t_jsb_pool)
      CALL mem%pools%Set_paths()

      ! add all variables in the sb_tissues_pool to the SB_ memory
      CALL mem%Add_var(mem%sb_main_pool%sb_pool, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%sb_main_pool%sb_formation, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes
      CALL mem%Add_var(mem%sb_main_pool%sb_loss, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes
      CALL mem%Add_var(mem%sb_main_pool%sb_transport, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes

      ! add all variables in the sb_elements_pool to the SB_ memory
      CALL mem%Add_var(mem%sb_main_pool%sb_delta_pools, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%sb_main_pool%sb_mycorrhiza_export, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes

      ! output to logfile all pools & their variables at model-init 
      ! (maybe use for debugging; its quite a long list)
      !IF (.NOT. ASSOCIATED(mem%parent)) THEN
      !  CALL mem%sb_main_pool%Print()
      !  ! CALL mem%sb_main_pool%sb_pool%Print()
      !END IF

      !=============================================================================================================================
      !> 2D & 3D variables
      !=============================================================================================================================

      CALL mem%Add_var('total_flux_carbon', mem%total_flux_carbon, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon', 'mol C m-2 s-1', 'carbon mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_nitrogen', mem%total_flux_nitrogen, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen', 'mol N m-2 s-1', 'nitrogen mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_phosphorus', mem%total_flux_phosphorus, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_phosphorus', 'mol P m-2 s-1', 'phosphorus mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_carbon13', mem%total_flux_carbon13, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon13', 'mol C13 m-2 s-1', 'carbon13 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_carbon14', mem%total_flux_carbon14, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon14', 'mol C14 m-2 s-1', 'carbon14 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_nitrogen15', mem%total_flux_nitrogen15, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen15', 'mol N15 m-2 s-1', 'nitrogen15 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('ph_sl', mem%ph_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('ph_sl', '', 'ph per soil layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue
    ENDIF

  END SUBROUTINE Init_sb_memory

END MODULE mo_sb_memory_class
